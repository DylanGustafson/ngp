//Non-generous prime searching program by Dylan G.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>         //only needed for memcpy
#include <gmp.h>            //GNU multiprecision arithmetic
#include <omp.h>            //OpenMP

#ifdef __APPLE__
  #include <sys/sysctl.h>   //MacOS-specific
#else
  #include <sys/sysinfo.h>  //Linux-specific
#endif

#define INIT_ALLOC 8  //Initial vector allocation: 6 prime factors (1 stdev above mean)
#define MID_ALLOC 12  //Mid allocation: 10 prime factors (3 stdev's above mean)
#define MAX_ALLOC 17  //Maximum possible (15) distinct prime factors for 64-bit value (ex: 693386350578511591)

typedef unsigned __int128 uint128_t;

//Global vars that remain constant across each thread
const char* primes_file = "primes.dat";
uint64_t chunk_size;
uint64_t n_classes;
uint32_t max_pfacs;
uint32_t* primes;
mpz_t* mpz_vars;
uint64_t* residue_counts;
uint64_t** vector_ptrs;

//Calculates a rough overestimate of n/ln(n)
uint64_t noln(uint64_t n) {
    //log2(n) approximated using clzll, 1/ln(2) approximated as 1.5
    n /= (63 - __builtin_clzll(n));
    return n + n / 2;
}

//Calculates a rough overestimate of sqrt(n)
uint32_t int_sqrt(uint64_t n) {
    uint32_t sn, mag;
    int clz = __builtin_clzll(n);

    if(!clz)
        //Linear approximation of sqrt(x) near 4294967296
        sn = 0x80000000 + (n >> 33);
    else {
        //Find closet power of 2 to sqrt(n), then do a Newton iteration
        mag = (64 - clz) / 2;
        sn = ((uint32_t)1 << (mag - 1)) + (uint32_t)(n >> (mag + 1));
    }
    return sn / 2 + n / sn / 2;   //One final iteration, avoiding overflow
}

//Reads in primes below 2^32 from primes.dat
void load_primes(uint64_t N_max) {
    uint32_t pfile_size, max_read, pmax;
    int mag;

    //fprintf(stderr, "Loading primes\n"); fflush(stderr);
    //Open primes.dat to get primes needed for sieve
    FILE* pfile = fopen(primes_file, "rb");
    if(pfile == NULL) {
        fprintf(stderr, "Error: could not find \"%s\" - run ./psieve first\n", primes_file);
        exit(1);
    }

    //First 4 bytes is total number of primes
    if(fread(&pfile_size, sizeof(uint32_t), 1, pfile) != 1) {
        fprintf(stderr, "Format error in \"%s\"\n", primes_file);
        exit(1);
    }

    //Allocate intial space for primes
    primes = (uint32_t*) malloc(pfile_size * sizeof(uint32_t));

    //Read pfile in 4KB blocks until primes are getting larger than sqrt(N_max)
    max_read = 0;
    pmax = int_sqrt(N_max);
    while(max_read < pfile_size) {
        max_read += fread(&primes[max_read], sizeof(uint32_t), 1024, pfile);

        //Break when final prime is greater or equal to sqrt(N_max)
        if(primes[max_read - 1] >= pmax)
            break;
    }
    fclose(pfile);

    //Find the actual index of first prime larger than pmax
    for(max_pfacs = max_read; max_pfacs >= max_read - 1024; --max_pfacs)
        if(primes[max_pfacs - 2] <= pmax)
            break;

    //Truncate primes array to free up memory for next steps
    primes = (uint32_t*) realloc(primes, max_pfacs * sizeof(uint32_t));
    //fprintf(stderr, "%u\n", primes[max_pfacs - 1]);
}

//Performs a sieve of eratosthenes on the chunk by setting vector pointers to null
void prime_sieve(const uint64_t chunk_start, uint64_t** vectors) {
    uint32_t i, p;
    uint64_t j, p2, mstart;

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

    //Copy in the tail end of the pattern up to the first multiple of 105
    uint64_t copy_start = (chunk_start + (pattern_size + 1) / 2) % pattern_size;
    if(__builtin_expect(copy_start != 0, 1)) {
        memcpy(&vectors[0], &pattern[copy_start], (pattern_size - copy_start) * sizeof(uint64_t*));
        copy_start = pattern_size - copy_start;
    }

    //Copy in the full pattern until it can't fit
    for(j = copy_start; j + pattern_size <= chunk_size; j += pattern_size)
        memcpy(&vectors[j], pattern, sizeof(pattern));

    //Copy in the beginning of the pattern up into the last bit of the chunk
    uint64_t tail_size = chunk_size - j;
    if(__builtin_expect(tail_size != 0, 1))
        memcpy(&vectors[j], pattern, tail_size * sizeof(uint64_t*));

    //Loop through each prime number, starting at 11
    for(i = 4; i < max_pfacs; ++i) {
        p = primes[i];
        p2 = (uint64_t) p * p / 2;

        //Get the offset of p^2 if it's in our range of N values
        // (Multiples of p below p^2 are already eliminated)
        mstart = p2 - chunk_start;

        //Otherwise just get the offset of the first N divisible by p
        if(p2 < chunk_start)
            mstart = p - 1 - (-mstart - 1) % p;

        //Set pointer to null in every position divisible by p
        if(p == 1) __builtin_unreachable();
        for(j = mstart; j < chunk_size; j += p)
            if(vectors[j])
                vectors[j] = 0;
    }
}

//For each N, fills up vectors with (N-1)/p values, where p is each prime factor of N-1
void gen_exp_vectors(const uint64_t chunk_start, uint64_t** vectors) {
    uint32_t i, p;
    uint64_t j, first_pos, N1_max, exp_val, newcap, divisor;
    uint128_t product;

    //Max 2n (ie N-1) value in chunk
    N1_max = 2 * (chunk_start + chunk_size);

    //Allocate vectors for each prime number in order to store (N-1)/p exponents
    for(j = 0; j < chunk_size; ++j) {
        if(!vectors[j])
            continue;

        vectors[j] = (uint64_t*) malloc(INIT_ALLOC * sizeof(uint64_t));
        vectors[j][0] = 3;              //Size
        vectors[j][1] = INIT_ALLOC;     //Capacity

        //First actual exponent value, which will correspond to the largest prime factor
        //of N-1. It's calculated by multiplying all prime factors less than sqrt(N-1).
        //Here we are using some bit trickery to take care of all powers of 2 right now
        vectors[j][2] = 2 * ((chunk_start + j) & -(chunk_start + j));
    }

    //Primary loop: Loop through the odd prime factors
    for(i = 1; i < max_pfacs; i++) {
        p = primes[i];

        //Get the offset of the first value divisible by p.
        first_pos = p - chunk_start % p;

        //The 2n/p (aka (N-1)/p) value just BEFORE the chunk
        //chunk_start / p has parentheses so that there is only one div instruction
        exp_val = 2 * (chunk_start / p);

        //Loop through chunk and push back exp_val into each vector
        if(p == 1) __builtin_unreachable();
        for(j = first_pos; j < chunk_size; j += p) {
            exp_val += 2;
            if(!vectors[j]) continue;

            //update fprod (first element of each vector)
            vectors[j][2] *= p;

            //To push back exp_val onto vector, first need to check if size == capacity
            if (vectors[j][0] == vectors[j][1]) {
                //Reallocate vector to increase size, and update capacity
                newcap = (vectors[j][1] == INIT_ALLOC) ? MID_ALLOC : MAX_ALLOC;
                vectors[j] = (uint64_t*) realloc(vectors[j], newcap * sizeof(uint64_t));
                vectors[j][1] = newcap;
            }

            //Place exponent value at the end of the vector and increment size
            vectors[j][vectors[j][0]] = exp_val;
            vectors[j][0] += 1;
        }

        //Skip next steps if p^2 already exceeds chunk's max N-1 value
        divisor = (uint64_t) p * p;
        if(divisor >= N1_max)
            continue;

        //Secondary Loop: Loop through powers of each prime p^e ('divisor').
        while(1) {
            //Get the offset of the first value divisible by p^e.
            first_pos = divisor - chunk_start % divisor;

            //Break if the first number divisble by p^e is outside of the current chunk.
            if(first_pos >= chunk_size)
                break;

            //Tertiary loop: Update fprods (first element of each vector)
            //do loop avoids an unnecessary conditional jump instruction
            j = first_pos;
            do {
                if(vectors[j])
                    vectors[j][2] *= p;

                //Saves space to help avoid cache misses
                if(divisor == 1) __builtin_unreachable();

                //j + divisor can't overflow unless chunk_size is on the order of 10^19
                //Example: try ./ngp 18401610824589482000 1000 1000 10
                j += divisor;
            } while(j < chunk_size);

            //Check if next p^e exceeds chunk's max N-1 value. Need this when there
            //is a very large p^e factor that divides an N-1 value in the chunk
            product = (uint128_t) divisor * p;
            if (product >= N1_max)
                break;

            //Product is stored in a dedicated uint128 before casting to uint64
            // to ensure that only one mul instruction is generated
            divisor = (uint64_t) product;
        }
    }
}

//Multiplication reduction subroutine for Montgomery modular multiplication
uint64_t mm_reduce(uint64_t a, uint64_t b, uint64_t p, uint64_t p_inv) {
    uint64_t result;

    uint128_t n = (uint128_t) a * b;

    //Calculate p * (n * p_inv % 2^64)
    uint128_t q = (uint128_t) p * (uint64_t)((uint64_t)n * p_inv);

    //arm64 specific instructions for calculating ((q + n) >> 64) % p
    #if defined(__aarch64__)
        asm volatile(
            "adds     x5, %[ql], %[nl]\n\t"     //128-bit addition of q and n, which stores
            "adcs     x6, %[qh], %[nh]\n\t"     // (q+n)>>64 into x6 and affects carry flag
            "sub      x5, x6, %[p]\n\t"         //Save ((q+n)>>64) - p into x5 in case it's needed
            "ccmp     x6, %[p], 2, cc\n\t"      //If q+n didnt carry, check if x6 > p, otherwise set carry
            "csel     %[res], x5, x6, hi"       //If x6 > p (or carry flag had been cleared), subtract p
            : [res] "=r" (result)
            : [ql] "r" ((uint64_t) q),
              [qh] "r" ((uint64_t) (q >> 64)),
              [nl] "r" ((uint64_t) n),
              [nh] "r" ((uint64_t) (n >> 64)),
              [p]  "r" (p)
            : "x5", "x6"                        //Since x0 - x3 are input args and x4 was used to calculate q
        );

        return result;

    //x86-64 specific instructions for calculating ((q + n) >> 64) % p
    #else
        asm volatile goto (
            "add     %[ql], %[nl]\n\t"          //128-bit addition of q and n, which
            "adc     %[qh], %[nh]\n\t"          // affects the carry flag if it overflows
            "mov     %[res], %[qh]\n\t"         //Store (q+n)>>64 into result var (%rax)
            "jc      %l[overflow]\n\t"          //If sum overflowed, jump to label
            "cmp     %[res], %[p]\n\t"          //Also need to jump there if
            "jnb     %l[overflow]"              // result is not below p
            : [res] "=r" (result)
            : [ql] "r" ((uint64_t) q),
              [qh] "r" ((uint64_t) (q >> 64)),
              [nl] "r" ((uint64_t) n),
              [nh] "r" ((uint64_t) (n >> 64)),
              [p]  "r" (p)
            : "cc"
            : overflow
        );
        return result;

    overflow:
        return result - p;

    #endif
}

//Finds the smallest positive primitive root of prime modulo using list of trial exponents.
uint64_t get_min_pr(const uint64_t p, const uint64_t* const exp_start, const uint64_t* const exp_end) {
    const uint64_t* ptr;
    uint64_t n, acc, res, root, root_mm;

    //2^128 % p
    uint64_t modp = 1 + (uint64_t)((uint128_t) -1 % p);

    //Calculate the inverse (modulo 2^64) of negative p
    uint64_t p_inv = (3 * -p) ^ 2;
    for(int i = 0; i < 4; ++i)
        p_inv *= (2 + p * p_inv);

    //Loop through root candidates, starting at 2
    root = 1;
    while(1) {
        //Incrementing at the start of the loop reduces the number of jump instructions
        ++root;

        //Convert root in preparation for Montgomery modular multiplication
        root_mm = mm_reduce(root, modp, p, p_inv);

        //Loop through trial exponents. Back to front is more efficient as it breaks sooner
        ptr = exp_start;
        while(1) {
            n = *ptr;

            // Perform modular exponentiation by squaring, with
            //  128-bit intermediate values (faster than GMP)
            res = 1;
            acc = root_mm;
            while(1) {
                if(n & 1)
                    res = mm_reduce(acc, res, p, p_inv);

                n >>= 1;
                if (n == 0) break;

                acc = mm_reduce(acc, acc, p, p_inv);
            }

            //Move on to next root candidate if result is ever 1
            if(res == 1) break;

            //Return the first root value that reaches the final exponent without breaking
            ++ptr;
            if(ptr == exp_end) return root;
        }
    }
}

//Thread entry point. Sieves primes, allcates exp vectors, calls the functions above
void find_mprs(const uint64_t chunk_start) {
    uint64_t j, modulo, root, chunk_j, vsize, fprod, residue;
    uint64_t *exp_start, *exp_end;
    int thread_num = omp_get_thread_num();

    //Pointer to this thread's array of exponent vector pointers
    uint64_t** vectors = vector_ptrs + thread_num * chunk_size;

    //Sieve out composite N (2n+1) values by making its vector pointer null
    //fprintf(stderr, "Sieving primes\n"); fflush(stderr);
    prime_sieve(chunk_start, vectors);

    //Fill vectors with (N-1)/p values for each prime number N
    //fprintf(stderr, "Factoring\n"); fflush(stderr);
    gen_exp_vectors(chunk_start, vectors);

    mpz_t* m1 = mpz_vars + 2 * thread_num;
    mpz_t* m2 = m1 + 1;
    uint64_t* res_counts_thr = residue_counts + thread_num * n_classes;

    //Get first primitive root of each prime N
    //fprintf(stderr, "Getting MPRs\n"); fflush(stderr);
    for(j = 0; j < chunk_size; ++j) {
        if(!vectors[j]) continue;

        //Calculate n, aka (N-1)/2
        chunk_j = chunk_start + j;

        //Get start and end pointers for exponents (the first two values are size and capacity)
        vsize = vectors[j][0];
        exp_start = vectors[j] + 2;
        exp_end = vectors[j] + vsize;

        //The first actual value in each vector ("fprod") is the product of all prime factors of N-1
        //(with exponents) below sqrt(N-1). If all prime factors of N-1 are less than sqrt(N-1), this
        //result in just N-1, so we dont need it. Since we do need (N-1)/2, we can write it there.
        if (*exp_start == chunk_j * 2) {
            *exp_start = chunk_j;
        }

        //If N-1 DOES have a prime factor larger than sqrt(N-1) (more likely), the result is now
        //(N-1)/p, where p is this largest prime factor. It's then ready to use as an exponent,
        //however we need to move it to the end and replace the first value with (N-1)/2, as it's
        //faster to check the exponents created from from smaller prime factors first.
        else {
            fprod = *exp_start;

            //If there is more capacity in the vector, simply move fprod to the
            //end, put (N-1)/2 in its original position, and increment the end pointer
            if(vectors[j][1] > vsize) {
                *exp_end = fprod;
                *exp_start = chunk_j;
                exp_end += 1;
            }

            //If there is no more capacity, move each value back one position
            //in the vector and put (N-1)/2 at the start and fprod at the end.
            //This overwrites the vector's capacity, but we don't need it anymore
            else {
                memmove(exp_start, exp_start + 1, (vsize - 3) * sizeof(uint64_t));
                *(exp_end - 1) = fprod;
                *(exp_start - 1) = chunk_j;
                exp_start -= 1;
            }
        }

        //Find first primitive root of N
        modulo = 2 * chunk_j + 1;
        root = get_min_pr(modulo, exp_start, exp_end);

        //Calculate residue: (root ^ (p - 1) % p^2) - 1
        mpz_set_ui(*m1, root);
        mpz_set_ui(*m2, modulo);
        mpz_mul_ui(*m2, *m2, modulo);
        mpz_powm_ui(*m1, *m1, modulo - 1, *m2);
        mpz_sub_ui(*m1, *m1, 1);

        //fprintf(stderr, "%19lu: %lu\n", modulo, root);
        //If result is 0, print out NGP
        if(!mpz_cmp_ui(*m1, 0))
            fprintf(stderr, "Non-generous prime found! Value = %lu (%lu)\n", modulo, root);

        //Calculate residue / p^2 * n_classes for 0 to 1 distribution
        mpz_mul_ui(*m1, *m1, n_classes);
        mpz_fdiv_q(*m1, *m1, *m2);
        residue = mpz_get_ui(*m1);
        res_counts_thr[residue] += 1;

        free(vectors[j]);
    }
}

//Main entry point. Reads in start, block_size and chunk_size from argv, then splits
//the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv) {
    char* endstr_ptr;
    uint64_t sys_memory, i, j;

    //Parse arguments, if any
    if(argc < 4) {
        puts("64-Bit Non-Generous Prime Searcher by Dylan G.\n");
        puts("Usage: ngp [c] start_value total_size chunk_size n_classes");
        puts("    c: Format output as csv, base 10");
        puts("    (default): tabular output, base 10");
        return 0;
    }

    //Read in numerical inputs
    uint64_t start = strtoull(argv[argc - 4], &endstr_ptr, 10) / 2;
    uint64_t block_size = strtoull(argv[argc - 3], &endstr_ptr, 10) / 2;
    chunk_size = strtoull(argv[argc - 2], &endstr_ptr, 10) / 2;
    n_classes = strtoull(argv[argc - 1], &endstr_ptr, 10);

    uint64_t N_min = 2 * start;
    uint64_t N_max = N_min + 2 * block_size;

    //Check if largest value overflows past 2^64
    if(N_max < N_min) {
        fprintf(stderr, "Error: maximum value exceeds 2^64-1\n");
        return 1;
    }

    //Check if chunk size is larger than total block size
    if(chunk_size > block_size) {
        fprintf(stderr, "Error: chunk size (%lu) must be smaller than total size (%lu)\n",
                chunk_size * 2, block_size * 2);
        return 1;
    }

    int max_threads = omp_get_max_threads();

    //Get RAM capacity in bytes using platform-specific system call
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

    //Estimate RAM usage using PNT. Shared: prime factors array. Thread: allocated vectors & ALL vector pointers.
    uint64_t primes_in_chunk = noln(N_min + 2 * chunk_size) - noln(N_min);
    uint64_t thread_mem = primes_in_chunk * MAX_ALLOC * sizeof(uint64_t) + chunk_size * sizeof(size_t) + n_classes * sizeof(uint64_t);

    uint64_t pfac_memory = noln(int_sqrt(N_max)) * sizeof(uint32_t);
    uint64_t req_memory = pfac_memory + thread_mem * max_threads;

    //Check RAM capacity against prediction
    if(req_memory > sys_memory) {
        fprintf(stderr, "Predicted RAM usage: %lu MiB (capacity: %lu MiB)\n", req_memory >> 20, sys_memory >> 20);
        fprintf(stderr, "Error: insufficient system memory\n");
        return 1;
    }

    //Load in primes up through sqrt(N_max) from primes.dat
    load_primes(N_max);

    //Allocate and initialize two mpz structs for each thread
    mpz_vars = (mpz_t*) malloc(2 * max_threads * sizeof(mpz_t));
    for(i = 0; i < 2 * max_threads; ++i)
        mpz_init(mpz_vars[i]);

    //Allocate residue arrays and vector pointers for each thread
    residue_counts = (uint64_t*) calloc(n_classes * max_threads, sizeof(uint64_t));
    vector_ptrs = (uint64_t**) malloc(chunk_size * max_threads * sizeof(uint64_t*));

    //Split chunks among threads using guided schedule since higher chunks contain fewer primes
    #pragma omp parallel for schedule(guided)
    for(i = 0; i < block_size / chunk_size; ++i)
        find_mprs(start + i * chunk_size);

    //Add up totals from each thread and format to stdout
    for(j = 0; j < n_classes; ++j) {
        for(i = 1; i < max_threads; ++i)
            residue_counts[j] += (residue_counts + i * n_classes)[j];

        if(argv[1][0] == 'c')
            printf("%lu,", residue_counts[j]);
        else
            printf("0.%03lu: %lu\n", j, residue_counts[j]);
    }

    //Garbage collect
    for(i = 0; i < 2 * max_threads; ++i)
        mpz_clear(mpz_vars[i]);

    free(mpz_vars);
    free(vector_ptrs);
    free(primes);

    return 0;
}
