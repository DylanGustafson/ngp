// Non-generous prime searching program by Dylan G.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>     // only needed for memcpy/memmove
#include <omp.h>        // OpenMP

// Number of bins used for q/N histogram
#define BIN_COUNT 1000

// Allocation size for vector of (N-1)/p exponents. 15 is the maximum
// number of prime factors for a 64-bit int (eg: 693386350578511591)
#define VECTOR_ALLOC 16

#define low64 0xFFFFFFFFFFFFFFFFULL

// Short aliases for fixed-width int types
typedef __uint128_t u128;
typedef uint64_t u64;
typedef uint32_t u32;

typedef struct{u64 x0; u64 x1; u64 x2; u64 x3;} u256;

// Globals (initialized before parallel section) which remain constant across threads
u64 block_start, block_size, chunk_size;
u64 vector_vals_size;
u32 max_pfacs;
u32* primes;
const char* primes_file = "primes.dat";

// Pointers to thread-partitioned arrays for safe parallel read/write
char* valid_shared;
u64** vectors_shared;
u64* vector_vals_shared;
u64* qcounts_shared;

// Converts numeric string to u64 type, while also handling integer scientific notation
u64 str_to_u64(char* str) {
    int overflow;
    int i = 0;
    u64 num = 0;
    const char err_msg[] = "Error: input larger than 2^64-1: %s\n";

    // Loop through digits until reaching 'e', 'E', or null char, or until overflow
    while ((str[i] | 0x20) != 'e') {
        if (str[i] == 0)
            return num;

        overflow = __builtin_mul_overflow(num, 10, &num);
        overflow |= __builtin_add_overflow(num, str[i] - '0', &num);
        if (overflow) {
            fprintf(stderr, err_msg, str);
            exit(1);
        }
        ++i;
    }

    // Read in 1 or 2 digit exponent
    int exp = str[i + 1] - '0';
    if (str[i + 2])
        exp = exp * 10 + str[i + 2] - '0';

    // Assume "e9" means "1e9"
    if (num == 0)
        num = 1;

    // Multiply by 10^exp and return
    for(i = 0; i < exp; ++i) {
        overflow = __builtin_mul_overflow(num, 10, &num);
        if (overflow) {
            fprintf(stderr, err_msg, str);
            exit(1);
        }
    }
    return num;
}

// Estimates the number of primes in [start, start + diff] using PNT
u64 estimate_nprimes(u64 start, u64 diff) {
    // log2(n) approximated using clzll, 1/ln(2) approximated as 1.5
    u64 result = diff + diff / 2;

    if (diff < start)
        result /= (63 - __builtin_clzll(start));
    else
        result /= (62 - __builtin_clzll(diff));

    //In extreme cases, just use 1000
    if (result < 1000)
        return 1000;

    return result;
}

// Calculates a rough overestimate of sqrt(n)
u32 int_sqrt(u64 n) {
    u32 sn, mag;
    int clz = __builtin_clzll(n);

    if (!clz)
        // Linear approximation of sqrt(x) near 4294967296
        sn = 0x80000000 + (n >> 33);
    else {
        // Find closet power of 2 to sqrt(n), then do a Newton iteration
        mag = (64 - clz) / 2;
        sn = ((u32)1 << (mag - 1)) + (u32)(n >> (mag + 1));
    }
    return sn / 2 + n / sn / 2;   // One final iteration, avoiding overflow
}

// Reads in primes below 2^32 from primes.dat
void load_primes(u64 N_max) {
    u32 pfile_size, max_read, pmax;

    // fprintf(stderr, "Loading primes\n"); fflush(stderr);
    // Open primes.dat to get primes needed for sieve
    FILE* pfile = fopen(primes_file, "rb");
    if (pfile == NULL) {
        fprintf(stderr, "Error: could not find \"%s\" - run ./psieve first\n", primes_file);
        exit(1);
    }

    // First 4 bytes is total number of primes
    if (fread(&pfile_size, sizeof(u32), 1, pfile) != 1) {
        fprintf(stderr, "Format error in \"%s\"\n", primes_file);
        exit(1);
    }

    // Allocate intial space for primes
    primes = (u32*) malloc(pfile_size * sizeof(u32));

    // Read pfile in 4KB blocks until primes are getting larger than sqrt(N_max)
    max_read = 0;
    pmax = int_sqrt(N_max);
    while (max_read < pfile_size) {
        max_read += fread(&primes[max_read], sizeof(u32), 1024, pfile);

        // Break when final prime is greater or equal to sqrt(N_max)
        if (primes[max_read - 1] >= pmax)
            break;
    }
    fclose(pfile);

    // Find the actual index of first prime larger than pmax
    for (max_pfacs = max_read; max_pfacs >= max_read - 1024; --max_pfacs)
        if (primes[max_pfacs - 2] <= pmax)
            break;

    // Truncate primes array to free up memory for next steps
    primes = (u32*) realloc(primes, max_pfacs * sizeof(u32));
    // fprintf(stderr, "%u\n", primes[max_pfacs - 1]);
}

// Performs a sieve of eratosthenes on the chunk by setting vector pointers to null
void prime_sieve(const u64 chunk_start, char* valid) {
    u32 i, p;
    u64 j, p2, mstart;

    const char pattern[] = {
        0,1,1,0,1,0,0,0,1,0,0,1,0,1,0,
        0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,
        0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
        0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,
        0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
        0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,
        0,0,1,0,1,0,0,1,0,0,0,1,0,1,1
    };

    const size_t pattern_size = sizeof(pattern);

    // Copy in the tail end of the pattern up to the first multiple of 105
    u64 copy_start = (chunk_start + (pattern_size + 1) / 2) % pattern_size;
    if (__builtin_expect(copy_start != 0, 1)) {
        memcpy(&valid[0], &pattern[copy_start], pattern_size - copy_start);
        copy_start = pattern_size - copy_start;
    }

    // Copy in the full pattern until it can't fit
    for (j = copy_start; j + pattern_size <= chunk_size; j += pattern_size)
        memcpy(&valid[j], pattern, sizeof(pattern));

    // Copy in the beginning of the pattern up into the last bit of the chunk
    u64 tail_size = chunk_size - j;
    if (__builtin_expect(tail_size != 0, 1))
        memcpy(&valid[j], pattern, tail_size);

    // Loop through each prime number, starting at 11
    for (i = 4; i < max_pfacs; ++i) {
        p = primes[i];
        p2 = (u64) p * p / 2;

        // Get the offset of p^2 if it's in our range of N values
        // (Multiples of p below p^2 are already eliminated)
        mstart = p2 - chunk_start;

        // Otherwise just get the offset of the first N divisible by p
        if (p2 < chunk_start)
            mstart = p - 1 - (-mstart - 1) % p;

        // Avoid unnecessary compiler optimization branch
        if (p == 1) __builtin_unreachable();

        // Set pointer to null in every position divisible by p
        for (j = mstart; j < chunk_size; j += p)
            if (valid[j])
                valid[j] = 0;
    }
}

// For each N, fills up vectors with (N-1)/p values, where p is each prime factor of N-1.
// Do loops are used throughout this function instead of for loops as they avoid the
// unnecessary cmp instruction generated at the top of each for loop, reducing binary size
void gen_exp_vectors(const u64 chunk_start, char* valid, u64** vectors, u64* vector_vals) {
    u32 i, p;
    u64 j, first_pos, N1_max, exp_val, divisor;
    u128 product;

    // Max 2n (ie N-1) value in chunk
    N1_max = 2 * (chunk_start + chunk_size);

    // "Allocate" space on vector_values for (N-1)/p exponents for each prime N
    u64* write_loc = vector_vals;
    j = 0;
    do {
        if (!valid[j])
            continue;

        vectors[j] = write_loc;
        write_loc += VECTOR_ALLOC;

        // Element 0 is vector size
        vectors[j][0] = 2;

        // First "actual" value will be HALF the (N-1)/p val corresponding to the largest prime
        // factor p. Thi is calculated via the product of all prime factors below sqrt(N-1),
        // with multiplicity. This initialization handles p=2 using a bitmasking trick
        vectors[j][1] = (chunk_start + j) & -(chunk_start + j);

    } while (++j != chunk_size);

    // If primes_in_chunk was an underestimate, exit and report required size
    if (__builtin_expect(write_loc - vector_vals > vector_vals_size, 0)) {
        fprintf(stderr, "Error: number of primes in chunk exceeds allocation %lu -> %lu ###\n",
                vector_vals_size / VECTOR_ALLOC,
                (u64) (write_loc - vector_vals) / VECTOR_ALLOC);
        exit(1);
    }

    // Primary loop: Loop through the odd prime factors
    i = 1;
    do {
        p = primes[i];

        // Get the offset of the first value divisible by p.
        first_pos = p - 1 - (chunk_start - 1) % p;

        // Move on to next prime if its first multiple is beyond the chunk
        if (first_pos >= chunk_size)
            continue;

        // The 2n/p (aka (N-1)/p) value just BEFORE the chunk
        // Must be written in this order to use the quotient value calculated in the
        // modulo operation above (%rax), thereby avoiding an extra div instruction
        exp_val = (chunk_start - 1) / p * 2;

        // Secondary Loop: Loop through chunk and push back exp_val into each vector
        j = first_pos;
        do {
            exp_val += 2;
            if (!valid[j]) continue;

            // Update fprod (first element of each vector)
            vectors[j][1] *= p;

            // Place exponent value at the end of the vector and increment size
            vectors[j][vectors[j][0]] = exp_val;
            vectors[j][0] += 1;
        } while ((j += p) < chunk_size);

        // Skip next steps if p^2 already exceeds chunk's max N-1 value
        divisor = (u64) p * p;
        if (divisor >= N1_max) continue;

        // Secondary Loop: Loop through powers of each prime p^e ('divisor').
        do {
            // Get the offset of the first value divisible by p^e.
            first_pos = divisor - 1 - (chunk_start - 1) % divisor;

            // Break if the first number divisble by p^e is outside of the current chunk.
            if (first_pos >= chunk_size) break;

            // Tertiary loop: Update fprods (first element of each vector)
            j = first_pos;
            do {
                if (valid[j])
                    vectors[j][1] *= p;

                // j + divisor can't overflow, even when chunk contains a multiple
                // of a large p^e value, unless chunk_size is on the order of 10^19
                j += divisor;
            } while (j < chunk_size);

            // Check if next p^e exceeds chunk's max N-1 value. Need 128 bits when
            // there is a large p^e factor that divides an N-1 value in the chunk
            product = (u128) divisor * p;

            // Product is stored in a dedicated uint128 before casting to uint64
            //  to ensure that only one mul instruction is generated
            divisor = (u64) product;

        } while (product < N1_max);
    } while (++i != max_pfacs);
}

// 64-bit REDC implementation for Montgomery modular multiplication
u64 redc64(u64 a, u64 b, u64 N, u64 N_inv) {
    u64 result;

    // m = N * (a * b * N_inv % 2^64), but need to save 128-bit product t
    u128 t = (u128) a * b;
    u128 m = (u128) N * (u64)((u64)t * N_inv);

    u64 th = (u64) (t >> 64);
    u64 mh = (u64) (m >> 64);

    // arm64 specific instructions for calculating ((m + 1) >> 64) % N
    #if defined(__aarch64__)
        asm volatile(
            "adds     x5, %[ml], %[tl]\n\t"     // 128-bit addition of m and t, which stores
            "adcs     x6, %[mh], %[th]\n\t"     // (m+t)>>64 into x6 and affects carry flag
            "sub      x5, x6, %[N]\n\t"         // Save ((m+t)>>64) - N into x5 in case it's needed
            "ccmp     x6, %[N], 2, cc\n\t"      // If m+t didnt carry, check if x6 > N, otherwise set carry
            "csel     %[res], x5, x6, hi"       // If x6 > N or carry flag had been cleared, subtract N

            : [res] "=r" (result)
            : [ml] "r" ((u64) m),
              [mh] "r" ((u64) (m >> 64)),
              [tl] "r" ((u64) t),
              [th] "r" ((u64) (t >> 64)),
              [N]  "r" (N)
            : "x5", "x6"                        // Since x0 - x3 are input args and x4 was used to calculate m
        );

        return result;

    // x86-64 specific instructions for calculating ((m - t) >> 64) % N
    #else
        asm volatile goto (
            "sub     %[mh], %[th]\n\t"
            "jnc     %l[noadd]"

            : [mh] "+r" (mh)
            : [th] "r" (th)
            : "cc"
            : noadd
        );

        return mh + N;
    noadd:
        return mh;

    /* Branchless version using cmov - slower!
        asm volatile(
            "sub     %[mh], %[th]\n\t"
            "lea     %[res], [%[mh]+%[N]]\n\t"
            "cmovnc  %[res], %[mh]"

            : [res] "=r" (result),
              [mh] "+r" (mh)
            : [th] "r" (th),
              [N]  "r" (N)
            : "cc"
        );
        return result;
    */

    #endif
}

// Finds the smallest positive primitive root of prime N using list of trial exponents.
u64 get_min_pr(const u64 N, const u64 N_inv, const u64 r128_N, const u64* const exp_start, const u64* const exp_end) {
    const u64* ptr;
    u64 exp, acc, res, root, start_mm;

    // Loop through root candidates, starting at 2
    root = 1;
    while(1) {
        // Increment root while skipping 4, 8, and 9. All perfect powers can be skipped,
        // but these are the only ones that save time. Incrementing root at the start of
        // the loop also saves time as it creates fewer total jump instructions.
        if ((root++ | 4) == 7)
            root += root / 4;

        // Convert root^2 in preparation for Montgomery modular multiplication
        start_mm = redc64(root * root, r128_N, N, N_inv);

        // Loop through (N-1)/p values, starting with smaller p values as it breaks sooner
        ptr = exp_start;
        while(1) {
            exp = *ptr;
            acc = start_mm;

            // Start with 1 unless exponent is odd, in which case need to start with root
            // This minor optimization is possible because (N-1)/p will never be 1
            res = (exp & 1) ? root : 1;
            exp >>= 1;

            // Perform modular exponentiation by squaring, with
            //  128-bit intermediate values (faster than GMP)
            while(1) {
                if (exp & 1)
                    res = redc64(acc, res, N, N_inv);

                exp >>= 1;
                if (exp == 0) break;

                acc = redc64(acc, acc, N, N_inv);
            }

            // Move on to next root candidate if result is ever 1
            if (res == 1) break;

            // Return the first root value that reaches the final exponent without breaking
            ++ptr;
            if (ptr == exp_end) return root;
        }
    }
}

// Calculates the 256-bit product of two unsigned 128-bit ints
u256 mult256(u128 a, u128 b) {
    u256 res;

    u64 al, ah, bl, bh;
    u64 carry3, carry4;
    u128 prd_00, prd_01, prd_10, prd_11, sum128;

    //split a and b into 64-bit halves
    al = (u64) a;
    ah = (u64) (a >> 64);
    bl = (u64) b;
    bh = (u64) (b >> 64);

    //Calculate 128-bit partial products
    prd_00 = (u128) al * bl;
    prd_01 = (u128) al * bh;
    prd_10 = (u128) ah * bl;
    prd_11 = (u128) ah * bh;

    res.x0 = (u64) prd_00;

    //128-bit sum of column 2, getting result & carry
    sum128 = (prd_00 >> 64) + (prd_01 & low64) + (prd_10 & low64);
    res.x1 = (u64) sum128;
    carry3 = (u64) (sum128 >> 64);

    //128-bit sum of column 3 (+ carry), getting result & carry
    sum128 = carry3 + (prd_01 >> 64) + (prd_10 >> 64) + (prd_11 & low64);
    res.x2 = (u64) sum128;
    carry4 = (u64) (sum128 >> 64);

    //64-bit sum of column 4 + carry, since it can't overflow
    res.x3 = carry4 + (u64)(prd_11 >> 64);

    return res;
}

//128-bit REDC implementation, no inline asm
/*
u128 redc128(u128 a, u128 b, u128 N2, u128 N2_inv) {
    u256 t = mult256(a, b);
    u128 t_lo = ((u128) t.x1 << 64) + t.x0;
    u128 t_hi = ((u128) t.x3 << 64) + t.x2;

    u256 m = mult256(N2, t_lo * N2_inv);
    u128 m_hi = ((u128) m.x3 << 64) + m.x2;

    u128 result = t_hi - m_hi;
    if (t_hi < m_hi)
        result += N2;

    return result;
}
*/

// 128-bit REDC implementation for Montgomery modular multiplication
u128 redc128(u128 a, u128 b, u128 N2, u128 N2_inv) {
    u256 t = mult256(a, b);

    // Construct 128-bit high & low words of product
    u128 t_lo = ((u128) t.x1 << 64) + t.x0;
    u128 t_hi = ((u128) t.x3 << 64) + t.x2;

    // Get the 128-bit high word of m = N2 * (t * N2_inv % 2^128)
    u256 m = mult256(N2, t_lo * N2_inv);
    u128 m_hi = ((u128) m.x3 << 64) + m.x2;

    // Decompose furter into 64-bit high-high & high-low words for inline asm below
    u64 t_hh = (u64) (t_hi >> 64);
    u64 t_hl = (u64)  t_hi;
    u64 m_hh = (u64) (m_hi >> 64);
    u64 m_hl = (u64)  m_hi;
    u64 N2_hi = (u64) (N2 >> 64);
    u64 N2_lo = (u64)  N2;

    // Subtract m_hi from t_hi, then add back N2 if it carried. This is
    // more efficient via inline asm as it avoids a needless 128-bit cmp
    asm volatile (
        "sub     %[thl], %[mhl]\n\t"
        "sbb     %[thh], %[mhh]\n\t"
        "jnc     1f\n\t"
        "add     %[thl], %[N2l]\n\t"
        "adc     %[thh], %[N2h]\n\t"
        "1:"
        : [thl] "+r" (t_hl), [thh] "+r" (t_hh)
        : [mhl] "r"  (m_hl), [mhh] "r"  (m_hh),
          [N2l] "r" (N2_lo), [N2h] "r" (N2_hi)
        : "cc"
    );

    // Rebuild the now-modified 128-bit t_hi
    return ((u128) t_hh << 64) + t_hl;
}

u64 calc_q(u64 N, const u64* exp_start, const u64* exp_end) {
    // Calculate modulo 2^128 multiplicative inverse of N
    u128 Ni_128 = ((u128) N * 3) ^ 2;
    u128 y = 1 - Ni_128 * N;

    // Five steps are required to guarantee 128 bits
    Ni_128 *= 1 + y;
    y *= y;
    Ni_128 *= 1 + y;
    y *= y;
    Ni_128 *= 1 + y;
    y *= y;
    Ni_128 *= 1 + y;
    y *= y;
    Ni_128 *= 1 + y;

    // Modular inverses needed for 64- and 128-bit montgomery form conversion
    u64 Ni_64 = (u64) Ni_128;
    u128 N2i_128 = Ni_128 * Ni_128;

    // Calculate 2^128 % N^2 as the first step towards the goal of 2^256 % N^2
    u128 N2 = (u128) N * N;
    u128 r128_N2 = 1 + (u128) -1 % N2;

    // Calculate 2^128 % N and (2^128 % N^2 - 2^128 % N) / N
    // using a single __udivmodti4 call
    u64 r128_N = r128_N2 % N;
    u64 quot = r128_N2 / N;

    // Since r128_N and quot are both less than N, we can choose whether
    // to multiply them or multiply their negatives modulo N. Using the
    // smaller value guarantees that multiplying by 2 wont overflow
    u128 prod1;
    if ((N - quot) < quot)
        prod1 = (u128) (N - r128_N) * ((N - quot) * 2);
    else
        prod1 = (u128) r128_N * (quot * 2);

    // Calculate prod2 = r128_N * quot * 2 * N % N2
    u64 mod1 = prod1 % N;
    u128 prod2 = (u128) mod1 * N;

    // Calculate the other product needed. Both products are less than N2.
    // so instead of a final modulo N2 we only need to potentially subtract
    // N2 once. However, since their sum can overflow past 2^128, its more
    // efficient to subtract them (need to use the negative mod N2) and just
    // add another N2 if the subtraction sets the carry flag
    u128 prod3 = N2 - (u128) r128_N * r128_N;

    // Decompose into 64-bit high & low words for inline assembly below
    u64 p2h = (u64) (prod2 >> 64);
    u64 p2l = (u64) prod2;
    u64 p3h = (u64) (prod3 >> 64);
    u64 p3l = (u64) prod3;
    u64 N2h = (u64) (N2 >> 64);
    u64 N2l = (u64)  N2;

    // Subtract prod3 from prod2, then add back N2 if it carried. This is
    // more efficient via inline asm as it avoids a needless 128-bit cmp
    asm volatile (
        "sub     %[p2l], %[p3l]\n\t"
        "sbb     %[p2h], %[p3h]\n\t"
        "jnc     1f\n\t"
        "add     %[p2l], %[N2l]\n\t"
        "adc     %[p2h], %[N2h]\n\t"
        "1:"
        : [p2l] "+r" (p2l), [p2h] "+r" (p2h)
        : [p3l]  "r" (p3l), [p3h]  "r" (p3h),
          [N2l]  "r" (N2l), [N2h]  "r" (N2h)
        : "cc"
    );

    // Construct 2^256 % N^2 from the 64-bit high & low words
    u128 r256_N2 = ((u128) p2h << 64) + p2l;

    // Get the minimum primitive root of N
    u64 root = get_min_pr(N, Ni_64, r128_N, exp_start, exp_end);

    // Set up 128-bit exponentiation. The exponent (N) is guaranteed
    // to be odd, so we can optimize past the first loop iteration
    u128 acc = redc128(root * root, r256_N2, N2, N2i_128);
    u128 res = root;
    u64 exp = N >> 1;

    // Use 128-bit montgomery exponentiation to calculate res = root^N % N^2
    while(1) {
        if(exp & 1)
            res = redc128(acc, res, N2, N2i_128);

        exp >>= 1;
        if(exp == 0) break;

        acc = redc128(acc, acc, N2, N2i_128);
    }

    // q = (root^N % N^2 - root) / N. Since this divides evenly, We
    // can use the previously calculated 128-bit modular inverse of N
    u64 q = (u64) ((res - root) * Ni_128);

    // Print "N (root q) factor_list" for debugging
    #ifdef verbose
    #pragma omp critical
    {
        fprintf(stderr, "%lu (%lu %lu)", N, root, q);
        for (int i = 0; i < exp_end - exp_start; ++i)
            fprintf(stderr, " %lu", (N - 1) / *(exp_start + i));
        fprintf(stderr, "\n");
    }
    #endif

    return q;
}

// Thread entry point. Sieves primes, allcates exp vectors, calls the functions above
void count_qvals(const u64 chunk_num) {
    u64 j, vsize, half_N, N, q_val;
    u64 *exp_start, *exp_end;
    u64 chunk_start = block_start + chunk_size * chunk_num;

    // Get pointers to thread-specific write locations
    int thread_num = omp_get_thread_num();
    char* valid = valid_shared + chunk_size * thread_num;
    u64** vectors = vectors_shared + chunk_size * thread_num;
    u64* vector_vals = vector_vals_shared + vector_vals_size * thread_num;
    u64* qcounts = qcounts_shared + BIN_COUNT * thread_num;

    #ifdef verbose
    #pragma omp critical
    {
        fprintf(stderr, "Thread %i (chunk %lu): Sieving\n", thread_num, chunk_num);
        fflush(stderr);
    }
    #endif

    // Sieve out composite N (2n+1) values by setting valid[] values to zero
    prime_sieve(chunk_start, valid);

    #ifdef verbose
    #pragma omp critical
    {
        fprintf(stderr, "Thread %i (chunk %lu): Factoring\n", thread_num, chunk_num);
        fflush(stderr);
    }
    #endif

    // Fill vectors with (N-1)/p values for each prime number N
    gen_exp_vectors(chunk_start, valid, vectors, vector_vals);

    #ifdef verbose
    #pragma omp critical
    {
        fprintf(stderr, "Thread %i (chunk %lu): Getting MPRs\n", thread_num, chunk_num);
        fflush(stderr);
    }
    #endif

    // Get first primitive root of each prime N and calculate q
    for (j = 0; j < chunk_size; ++j) {
        if (!valid[j]) continue;

        // Get start and end pointers for exponents
        vsize = vectors[j][0];
        exp_start = vectors[j] + 1;
        exp_end = vectors[j] + vsize;

        // Calculate (N-1)/2
        half_N = chunk_start + j;

        // The first actual value in each vector ("fprod") is HALF the product of all prime factors
        // of N-1 (with multiplicity) below sqrt(N-1). If ALL prime factors of N-1 are less than
        // sqrt(N-1), this result is just (N-1)/2, which can now be used as the first exponent (p=2).
        // If N-1 DOES have a prime factor larger than sqrt(N-1) (more likely), the result is now
        // (N-1)/(p*2), where p is this largest prime factor. We need to double this value and move
        // it to the end of the vector, and replace the first value with (N-1)/2. When looking for the
        // MPR, larger (N-1)/p values fail sooner (r^e % p = 1), so they should be sorted descending.
        if (*exp_start != half_N) {
            *exp_end = *exp_start * 2;
            exp_end += 1;
            *exp_start = half_N;
        }

        N = 2 * half_N + 1;
        q_val = calc_q(N, exp_start, exp_end);

        // Print numbers of interest to stderr
        if (__builtin_expect(q_val == 0, 0))
            #pragma omp critical
            fprintf(stderr, "Non-generous prime found! Value = %lu\n", N);
        else if (__builtin_expect(q_val == 1, 0))
            #pragma omp critical
            fprintf(stderr, "q = 1. Value = %lu\n", N);

        // Increment counter in q/N bin
        qcounts[(u128)q_val * BIN_COUNT / N] += 1;
    }

    #ifdef verbose
    #pragma omp critical
    {
        fprintf(stderr, "Thread %i (chunk %lu): Done\n", thread_num, chunk_num);
        fflush(stderr);
    }
    #endif
}

// Main entry point. Reads in start, block_size and chunk_size from argv, then splits
// the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv) {
    u64 req_mem, primes_in_chunk, i, j;

    // Parse arguments, if any
    if (argc < 4) {
        puts("64-Bit Non-Generous Prime Searcher by Dylan G.\n");
        puts("Usage: ngp [m|c] [p=nprimes] start_value total_size chunk_size");
        puts("    m: Only return predicted memory usage for given args");
        puts("    c: Format output as csv - defaults to tabular output");
        puts("    p: Specify number of primes per chunk (for malloc) - defaults to PNT estimate");
        return 0;
    }

    // Read in numerical inputs
    block_start = str_to_u64(argv[argc - 3]) / 2;
    block_size = str_to_u64(argv[argc - 2]) / 2;
    chunk_size = str_to_u64(argv[argc - 1]) / 2;

    u64 N_min = 2 * block_start;
    u64 N_max = N_min + 2 * block_size;

    // Check if largest value overflows past 2^64
    if (N_max < N_min) {
        fprintf(stderr, "Error: maximum value exceeds 2^64-1\n");
        return 1;
    }

    // Check if chunk size is larger than total block size
    if (chunk_size > block_size) {
        fprintf(stderr, "Error: chunk size (%lu) must be smaller than total size (%lu)\n",
                chunk_size * 2, block_size * 2);
        return 1;
    }

    // Handle starting values below 10
    int add_2357 = 0;
    if (block_start < 5) {
        block_start = 5;
        block_size -= 5;
        add_2357 = 1;

        // For 1 trillion, split into 974205 chunks
        if (block_size == 1e12 / 2 - 5)
            chunk_size = 513239;

        // For 1 billion, split into 10001 chunks
        else if (block_size == 1e9 / 2 - 5)
            chunk_size = 49995;

        //For 10 million or less, just use a single thread
        else if (block_size < 1e7 / 2)
            chunk_size = block_size;

        // Otherwise exit
        else {
            fprintf(stderr, "start_value must be >=10, OR use total_size of 1e9, 1e12, or <=1e7\n");
            return 1;
        }

        fprintf(stderr, "Changed args to %lu %lu %lu (p < 10 data will be added later)\n",
                block_start * 2, block_size * 2, chunk_size * 2);
    }

    // Check if there will be an integer number of chunks
    u64 num_chunks = block_size / chunk_size;
    if (block_size % chunk_size != 0) {
        fprintf(stderr, "Error: chunk_size must evenly divide total_size\n");
        return 1;
    }

    // If number of primes in chunk was not given, estimate using PNT
    if (argv[argc - 4][0] == 'p')
        primes_in_chunk = str_to_u64(&argv[argc - 4][2]);
    else
        primes_in_chunk = estimate_nprimes(N_min, chunk_size * 2);

    vector_vals_size = primes_in_chunk * VECTOR_ALLOC;

    // Predict memory usage
    int max_threads = omp_get_max_threads();
    req_mem = 8e8 + max_threads * (
        BIN_COUNT * sizeof(u64) +
        chunk_size * (sizeof(u64*) + 1) +
        vector_vals_size * sizeof(u64));

    if (argv[1][0] == 'm') {
        fprintf(stderr, "Predicted usage: %lu MiB\n", req_mem >> 20);
        return 0;
    }

    // Load in prime factors into global primes[] (read-only for threads)
    load_primes(N_max);

    // Allocate partitioned read/write arrays for each thread
    valid_shared = (char*) aligned_alloc(64, max_threads * chunk_size);
    vectors_shared = (u64**) malloc(max_threads * chunk_size * sizeof(u64*));
    vector_vals_shared = (u64*) aligned_alloc(64, max_threads * vector_vals_size * sizeof(u64));
    qcounts_shared = (u64*) calloc(max_threads * BIN_COUNT, sizeof(u64));

    // Exit if any allocation(s) failed
    if (!qcounts_shared || !valid_shared || !vectors_shared || !vector_vals_shared) {
        fprintf(stderr, "Error: couldn't allocate necessary memory (%lu GiB)\n", req_mem >> 30);
        return 1;
    }

    // Split chunks among threads using guided schedule since higher chunks contain fewer primes
    #pragma omp parallel for schedule(guided)
    for (i = 0; i < num_chunks; ++i)
        count_qvals(i);

    // If 0 - 10 was skipped, add in hard-coded q/N values for 2, 3, 5, and 7
    if (add_2357) {
        qcounts_shared[0]++;
        qcounts_shared[2 * BIN_COUNT / 3]++;
        qcounts_shared[1 * BIN_COUNT / 5]++;
        qcounts_shared[4 * BIN_COUNT / 7]++;
    }

    // Add up totals from each thread and format to stdout
    for (j = 0; j < BIN_COUNT; ++j) {
        for (i = 1; i < max_threads; ++i)
            qcounts_shared[j] += (qcounts_shared + i * BIN_COUNT)[j];

        if (argv[1][0] == 'c')
            printf("%lu,", qcounts_shared[j]);
        else
            printf("0.%03lu: %lu\n", j, qcounts_shared[j]);
    }

    // Garbage collect
    free(qcounts_shared);
    free(vector_vals_shared);
    free(vectors_shared);
    free(valid_shared);
    free(primes);

    return 0;
}
