// Non-generous prime searching program by Dylan G.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>        // only needed for memcpy
#include <gmp.h>           // GNU multiprecision arithmetic
#include <omp.h>           // OpenMP

#ifdef __APPLE__
  #include <sys/sysctl.h>  // MacOS-specific
#else
  #include <sys/sysinfo.h> // Linux-specific
#endif

#define BIN_COUNT 1000     // Number of bins used for q/N histogram
#define INIT_ALLOC 8       // Initial vector allocation: 6 prime factors (1 stdev above mean)
#define MID_ALLOC 12       // Mid allocation: 10 prime factors (3 stdev's above mean)
#define MAX_ALLOC 17       // Maximum possible (15) distinct prime factors for 64-bit int (eg: 693386350578511591)

typedef __uint128_t u128;
typedef uint64_t u64;
typedef uint32_t u32;

// Global vars that remain constant across each thread
const char* primes_file = "primes.dat";
u64 chunk_size;
u32 max_pfacs;
u32* primes;
mpz_t* mpz_vars;
u64* residue_counts;
u64** vector_ptrs;

// Converts numeric string to u64 type, while also handling integer scientific notation
u64 str_to_u64(char* str) {
    int i, overflow;
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

// Calculates a rough overestimate of n/ln(n)
u64 noln(u64 n) {
    // log2(n) approximated using clzll, 1/ln(2) approximated as 1.5
    n /= (63 - __builtin_clzll(n));
    return n + n / 2;
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
void prime_sieve(const u64 chunk_start, u64** vectors) {
    u32 i, p;
    u64 j, p2, mstart;

    const size_t pattern[] = {
        0,1,1,0,1,0,0,0,1,0,0,1,0,1,0,
        0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,
        0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
        0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,
        0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
        0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,
        0,0,1,0,1,0,0,1,0,0,0,1,0,1,1,
    };

    const size_t pattern_size = sizeof(pattern) / sizeof(pattern[0]);

    // Copy in the tail end of the pattern up to the first multiple of 105
    u64 copy_start = (chunk_start + (pattern_size + 1) / 2) % pattern_size;
    if (__builtin_expect(copy_start != 0, 1)) {
        memcpy(&vectors[0], &pattern[copy_start], (pattern_size - copy_start) * sizeof(u64*));
        copy_start = pattern_size - copy_start;
    }

    // Copy in the full pattern until it can't fit
    for (j = copy_start; j + pattern_size <= chunk_size; j += pattern_size)
        memcpy(&vectors[j], pattern, sizeof(pattern));

    // Copy in the beginning of the pattern up into the last bit of the chunk
    u64 tail_size = chunk_size - j;
    if (__builtin_expect(tail_size != 0, 1))
        memcpy(&vectors[j], pattern, tail_size * sizeof(u64*));

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
            if (vectors[j])
                vectors[j] = 0;
    }
}

// For each N, fills up vectors with (N-1)/p values, where p is each prime factor of N-1.
// Do loops are used throughout this function instead of for loops as they avoid the
// unnecessary cmp instruction generated at the top of each for loop, reducing binary size
void gen_exp_vectors(const u64 chunk_start, u64** vectors) {
    u32 i, p;
    u64 j, first_pos, N1_max, exp_val, newcap, divisor;
    u128 product;

    // Max 2n (ie N-1) value in chunk
    N1_max = 2 * (chunk_start + chunk_size);

    // Allocate vectors for each prime number in order to store (N-1)/p exponents
    j = 0;
    do {
        if (!vectors[j])
            continue;

        vectors[j] = (u64*) malloc(INIT_ALLOC * sizeof(u64));
        vectors[j][0] = 3;           // Size
        vectors[j][1] = INIT_ALLOC;  // Capacity

        // First "actual" exponent value, which will correspond to the largest prime
        // factor of N-1. It's calculated by multiplying all prime factors (with
        // multiplicity) less than sqrt(N-1). This initialization expression uses
        // bit manipulation to handle p=2 by calculating the largest 2^(e-1) factor
        vectors[j][2] = (chunk_start + j) & -(chunk_start + j);

    } while (++j != chunk_size);

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
        // Must be written in this order so it can use the %rax
        // value (quotient) calculated in the modulo operation above
        exp_val = (chunk_start - 1) / p * 2;

        // Secondary Loop: Loop through chunk and push back exp_val into each vector
        j = first_pos;
        do {
            exp_val += 2;
            if (!vectors[j]) continue;

            // Update fprod (first element of each vector)
            vectors[j][2] *= p;

            // If vector's size == capacity, reallocate with larger capacity
            if (vectors[j][0] == vectors[j][1]) {
                newcap = (vectors[j][1] == INIT_ALLOC) ? MID_ALLOC : MAX_ALLOC;
                vectors[j] = (u64*) realloc(vectors[j], newcap * sizeof(u64));
                vectors[j][1] = newcap;
            }

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
                if (vectors[j])
                    vectors[j][2] *= p;

                // j + divisor can't overflow, even when chunk contains a multiple
                // of a large p^e value, unless chunk_size is on the order of 10^19
                // Example: try ./ngp 18401610824589482000 1000 1000
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

// 64-bit REDC subroutine for Montgomery modular multiplication
u64 redc64(u64 a, u64 b, u64 N, u64 N_inv) {
    u64 result;

    // m = N * (a * b * N_inv % 2^64), but need to save 128-bit product t
    u128 t = (u128) a * b;
    u128 m = (u128) N * (u64)((u64)t * N_inv);

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

    // x86-64 specific instructions for calculating ((m + t) >> 64) % N
    #else
        asm volatile goto (
            "add     %[ml], %[tl]\n\t"          //  128-bit addition of m and t, which
            "adc     %[mh], %[th]\n\t"          // affects the carry flag if it overflows
            "mov     %[res], %[mh]\n\t"         // Store (m+t)>>64 into result var (%rax)
            "jc      %l[over]\n\t"              // If sum overflowed, jump to label
            "cmp     %[res], %[N]\n\t"          // Also need to jump there if
            "jnb     %l[over]"                  //  result is not below p
            : [res] "=r" (result)
            : [ml] "r" ((u64) m),
              [mh] "r" ((u64) (m >> 64)),
              [tl] "r" ((u64) t),
              [th] "r" ((u64) (t >> 64)),
              [N]  "r" (N)
            : "cc"
            : over
        );
        return result;

    over:
        return result - N;

    #endif
}

// Finds the smallest positive primitive root of prime N using list of trial exponents.
u64 get_min_pr(const u64 N, const u64* const exp_start, const u64* const exp_end) {
    const u64* ptr;
    u64 exp, acc, res, root, start_mm;

    // 2^128 % N
    u64 R2modN = 1 + (u64) ((u128) -1 % N);

    // Calculate the inverse (modulo 2^64) of negative N
    u64 N_inv = (3 * -N) ^ 2;
    for (int i = 0; i < 4; ++i)
        N_inv *= (2 + N * N_inv);

    // Loop through root candidates, starting at 2
    root = 1;
    while(1) {
        // Increment root while skipping 4, 8, and 9. All perfect powers can be skipped,
        // but these are the only ones that save time. Incrementing root at the start of
        // the loop also saves time as it creates fewer total jump instructions.
        if ((root++ | 4) == 7)
            root += root / 4;

        // Convert root^2 in preparation for Montgomery modular multiplication
        start_mm = redc64(root * root, R2modN, N, N_inv);

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

// Thread entry point. Sieves primes, allcates exp vectors, calls the functions above
void find_mprs(const u64 chunk_start) {
    u64 j, N, root, half_N, vsize, fprod, q, index;
    u64 *exp_start, *exp_end;

    // Get pointers to thread-specific write locations
    int thread_num = omp_get_thread_num();
    u64** vectors = vector_ptrs + chunk_size * thread_num;
    u64* res_counts_thr = residue_counts + BIN_COUNT * thread_num;
    mpz_t* m1 = mpz_vars + 2 * thread_num;
    mpz_t* m2 = m1 + 1;

    // Sieve out composite N (2n+1) values by making its vector pointer null
    //fprintf(stderr, "Sieving primes\n"); fflush(stderr);
    prime_sieve(chunk_start, vectors);

    // Fill vectors with (N-1)/p values for each prime number N
    //fprintf(stderr, "Factoring\n"); fflush(stderr);
    gen_exp_vectors(chunk_start, vectors);

    // Get first primitive root of each prime N and calculate q
    //fprintf(stderr, "Getting MPRs\n"); fflush(stderr);
    for (j = 0; j < chunk_size; ++j) {
        if (!vectors[j]) continue;

        // Calculate (N-1)/2
        half_N = chunk_start + j;

        // Get start and end pointers for exponents (the first two values are size and capacity)
        vsize = vectors[j][0];
        exp_start = vectors[j] + 2;
        exp_end = vectors[j] + vsize;

        // The first actual value in each vector ("fprod") is HALF the product of all prime factors
        // of N-1 (with multiplicity) below sqrt(N-1). If all prime factors of N-1 are less than
        // sqrt(N-1), this result is (N-1)/2, so we can just use it as the first exponent (p=2).
        // If N-1 DOES have a prime factor larger than sqrt(N-1) (more likely), the result is now
        // (N-1)/(p*2), where p is this largest prime factor. We need to double this value and move
        // it to the end of the vector, and replace the first value with (N-1)/2, as it's faster to
        // check the exponents created from from smaller prime factors first.
        if (*exp_start != half_N) {
            fprod = *exp_start * 2;

            // If the vector has capacity, move fprod to the end and increment the end pointer
            if (vectors[j][1] > vsize) {
                *exp_end = fprod;
                exp_end += 1;
            }

            // If there is no more capacity, move each value back one position
            // in the vector and put (N-1)/2 at the start and fprod at the end.
            // This overwrites the vector's capacity, but we don't need it anymore
            else {
                memmove(exp_start, exp_start + 1, (vsize - 3) * sizeof(u64));
                *(exp_end - 1) = fprod;
                exp_start -= 1;
            }

            *exp_start = half_N;
        }

        // Find first primitive root of N
        N = 2 * half_N + 1;
        root = get_min_pr(N, exp_start, exp_end);

        // Calculate q = (root ^ N % N^2 - root) / N
        mpz_set_ui(*m1, root);
        mpz_set_ui(*m2, N);
        mpz_mul_ui(*m2, *m2, N);
        mpz_powm_ui(*m1, *m1, N, *m2);
        mpz_sub_ui(*m1, *m1, root);
        mpz_divexact_ui(*m1, *m1, N);
        q = mpz_get_ui(*m1);

        //  Print numbers of interest to stderr
        if (__builtin_expect(q == 0, 0))
            fprintf(stderr, "Non-generous prime found! Value = %lu (%lu)\n", N, root);
        else if (__builtin_expect(q == 1, 0))
            fprintf(stderr, "q = 1. Value = %lu (%lu)\n", N, root);

        // Verbose output: Print N, root, q, and prime factors of N-1 to stderr
        #ifdef verbose
        fprintf(stderr, "%lu (%lu %lu)", N, root, q);
        for (int i = 0; i < exp_end - exp_start; ++i)
            fprintf(stderr, " %lu", (N - 1) / *(exp_start + i));
        fprintf(stderr, "\n");
        #endif

        // Increment counter in q/N bin
        res_counts_thr[(u128)q * BIN_COUNT / N] += 1;
        free(vectors[j]);
    }
}

// Main entry point. Reads in start, block_size and chunk_size from argv, then splits
// the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv) {
    u64 sys_memory, i, j;
    int add_2357 = 0;

    // Parse arguments, if any
    if (argc < 4) {
        puts("64-Bit Non-Generous Prime Searcher by Dylan G.\n");
        puts("Usage: ngp [c] start_value total_size chunk_size");
        puts("    c: Format output as csv, base 10");
        puts("    (default): tabular output, base 10");
        return 0;
    }

    // Read in numerical inputs
    u64 block_start = str_to_u64(argv[argc - 3]) / 2;
    u64 block_size = str_to_u64(argv[argc - 2]) / 2;
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
        fprintf(stderr, "Error: chunk size (%lu) must be smaller than total size (%lu)\n", chunk_size * 2, block_size * 2);
        return 1;
    }

    // Handle starting values below 10
    if (block_start < 5) {
        block_start = 5;
        block_size -= 5;
        add_2357 = 1;

        // For 1 billion, split into 10001 chunks
        if (block_size == 499999995)
            chunk_size = 49995;

        //For 10 million or less, just use a single thread
        else if (block_size < 5000000)
            chunk_size = block_size;

        // Otherwise exit
        else {
            fprintf(stderr, "start_value must be at least 10, OR use total_size of 1e9 or <1e7\n");
            return 1;
        }

        fprintf(stderr, "Changed args to %lu %lu %lu\n", block_start * 2, block_size * 2, chunk_size * 2);
    }

    // Get RAM capacity in bytes using platform-specific system call
    #ifdef __APPLE__
        int mib[2] = {CTL_HW, HW_MEMSIZE};
        size_t len_memory = sizeof(sys_memory);
        
        if (sysctl(mib, 2, &sys_memory, &len_memory, NULL, 0) == -1) {
            perror("sysctl");
            return 1;
        }
    #else
        struct sysinfo info;
        if (sysinfo(&info) != 0) {
            perror("sysinfo");
            return 1;
        }
        sys_memory = info.totalram;
    #endif

    // Estimate RAM usage using PNT. Shared: prime factors array. Thread: allocated vectors & ALL vector pointers.
    int max_threads = omp_get_max_threads();
    u64 primes_in_chunk = noln(N_min + 2 * chunk_size) - noln(N_min);
    u64 thread_mem = primes_in_chunk * MAX_ALLOC * sizeof(u64) + chunk_size * sizeof(size_t) + BIN_COUNT * sizeof(u64);
    u64 pfac_memory = noln(int_sqrt(N_max)) * sizeof(u32);
    u64 req_memory = pfac_memory + thread_mem * max_threads;

    // Check RAM capacity against prediction
    if (req_memory > sys_memory) {
        fprintf(stderr, "Error: insufficient system memory\n");
        fprintf(stderr, "Predicted RAM usage: %lu MiB (capacity: %lu MiB)\n", req_memory >> 20, sys_memory >> 20);
        return 1;
    }

    // Load in primes up through sqrt(N_max) from primes.dat
    load_primes(N_max);

    // Allocate and initialize two mpz structs for each thread
    mpz_vars = (mpz_t*) malloc(2 * max_threads * sizeof(mpz_t));
    for (i = 0; i < 2 * max_threads; ++i)
        mpz_init(mpz_vars[i]);

    // Allocate residue arrays and vector pointers for each thread
    residue_counts = (u64*) calloc(BIN_COUNT * max_threads, sizeof(u64));
    vector_ptrs = (u64**) malloc(chunk_size * max_threads * sizeof(u64*));

    // Split chunks among threads using guided schedule since higher chunks contain fewer primes
    #pragma omp parallel for schedule(guided)
    for (i = 0; i < block_size / chunk_size; ++i)
        find_mprs(block_start + i * chunk_size);

    // If 0 - 10 was skipped, add in hard-coded q/N values for 2, 3, 5, and 7
    if (add_2357) {
        residue_counts[0]++;
        residue_counts[2 * BIN_COUNT / 3]++;
        residue_counts[1 * BIN_COUNT / 5]++;
        residue_counts[4 * BIN_COUNT / 7]++;
    }

    // Add up totals from each thread and format to stdout
    for (j = 0; j < BIN_COUNT; ++j) {
        for (i = 1; i < max_threads; ++i)
            residue_counts[j] += (residue_counts + i * BIN_COUNT)[j];

        if (argv[1][0] == 'c')
            printf("%lu,", residue_counts[j]);
        else
            printf("0.%03lu: %lu\n", j, residue_counts[j]);
    }

    // Garbage collect
    for (i = 0; i < 2 * max_threads; ++i)
        mpz_clear(mpz_vars[i]);

    free(mpz_vars);
    free(vector_ptrs);
    free(primes);

    return 0;
}
