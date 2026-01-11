//Non-generous prime searching program by Dylan G.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>         //only needed for memmove
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

//Output structure for each prime
typedef struct {
    uint64_t prime;
    uint64_t root;
    uint64_t residue;
} archive;

//IO file names
const char* primes_file = "primes.dat";

//Global vars that remain constant across each thread
uint64_t chunk_size;
uint32_t* primes;
uint32_t max_pfacs;
uint32_t output_format;

//256-bit values private to each thread
mpz_t g_a, g_b;
#pragma omp threadprivate(g_a)
#pragma omp threadprivate(g_b)

//Overstimates n/ln(n) using simple integer math
uint64_t noln(uint64_t n) {
    //log2(n) approximated using clzll, 1/ln(2) approximated as 1.5
    n /= (64 - __builtin_clzll(n));
    return n + n / 2;
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

    //Calculate pmax as a rough overestimate of sqrt(N_max)
    if(!__builtin_clzll(N_max))
        //Linear approximation of sqrt(x) near 4294967296
        pmax = 0x80000000 + (N_max >> 33);
    else {
        //Find closet power of 2 to sqrt(N_max), then do a Newton iteration
        mag = (64 - __builtin_clzll(N_max)) / 2;
        pmax = ((uint32_t)1 << (mag - 1)) + (uint32_t)(N_max >> (mag + 1));
    }
    pmax = pmax / 2 + N_max / pmax / 2;   //One final iteration, avoiding overflow

    //Read pfile in 4KB blocks until primes are getting larger than pmax
    max_read = 0;
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
    uint64_t j, p2, N_min, mstart;

    //First 2n+1 value for calculating offsets
    N_min = 2 * chunk_start + 1;

    //Loop through each prime number except 2 (since 2n+1 is always odd)
    for(i = 1; i < max_pfacs; ++i) {
        p = primes[i];
        p2 = (uint64_t) p * p;

        //Get the offset of p^2 if it's in our range of N values
        //(Multiples of p below p^2 are already eliminated)
        if (p2 > N_min)
            mstart = (p2 - N_min) / 2;

        //Otherwise just get the offset of the first N divisible by p
        else
            mstart = p - 1 - ((N_min - p) / 2 - 1) % p;

        //Set pointer to null in every position divisible by p
        for(j = mstart; j < chunk_size; j += p)
            if(vectors[j])
                vectors[j] = 0;
    }
}

//For each N, fills up vectors with (N-1)/p values, where p is each prime factor of N-1
void gen_exps(const uint64_t chunk_start, uint64_t** vectors) {
    uint32_t i, p;
    uint64_t j, first_pos, N1_max, exp_val, newcap, divisor;
    uint128_t product;

    N1_max = 2 * (chunk_start + chunk_size);

    //Primary loop: Loop through prime factors (p = 2 was already accounted for)
    for(i = 1; i < max_pfacs; i++) {
        p = primes[i];

        //Get the offset of the first value divisible by p.
        first_pos = p - chunk_start % p;

        //First 2n/p (aka (N-1)/p) value, minus 2 to fit with for loop
        exp_val = 2 * (chunk_start + first_pos) / p - 2;

        //Loop through chunk and push back exp_val into each vector
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
            j = first_pos;
            do {
                if(vectors[j])
                    vectors[j][2] *= p;

                //This sum cant overflow unless chunk_size is on the order of 10^19
                //Proof: try ngp 18401610824589482000 1000 1000
                j += divisor;
            } while(j < chunk_size);

            //Check if next p^e exceeds chunk's max N-1 value. Need this when there
            //is a very large p^e factor that divides an N-1 value in the chunk
            product = (uint128_t) divisor * p;
            if (product >= N1_max)
                break;

            //Product is stored in a dedicated uint128 before casting to uint64
            //to ensure that only one mul instruction is generated
            divisor = (uint64_t) product;
        }
    }
}

//Calculates the modular inverse (base 2^64) of a
uint64_t mod_inv64(uint64_t a) {
    uint64_t r1, r2, s1, s2, quot, swap;

    r2 = a;
    s1 = 0;
    s2 = 1;

    //Used instead of 2^64, works the same!
    quot = 0xFFFFFFFFFFFFFFFF / r2;
    swap = -quot * r2;

    //Do loop is used to avoid an unecessary cmp instruction
    do {
        //Extended Euclidean algorithm
        r1 = r2;
        r2 = swap;

        swap = s1 - quot * s2;
        s1 = s2;
        s2 = swap;

        quot = r1 / r2;
        swap = r1 - quot * r2;
    } while(swap);

    return s2;
}

//Multiplication reduction subroutine for Montgomery modular multiplication
uint64_t mm_reduce(uint128_t n, uint64_t p, uint64_t p_i) {
    uint64_t result;

    //Calculate p * (n * p_i % 2^64)
    uint128_t q = (uint128_t) p * (uint64_t)(n * p_i);

    //arm64-specific instructions for efficiently calculating ((q + n) >> 64) % p
    #if defined(__aarch64__)
        asm volatile(
            "adds     x5, %[ql], %[nl]\n\t"     // 128-bit addition of q and n, which stores
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

    //x86-64 specific instructions for ((q + n) >> 64) % p: just need to save carry flag into %eax
    #else
        //Adding q and n. This 128-bit sum can overflow!
        uint128_t sum = q + n;

        //Need to save the state of the carry flag resulting from the sum above.
        //The sbb instruction sets eax to 0 if C is not set, but -1 if C is set.
        //By storing this result in chk_carry, we can include this as a necessary
        //condition to subtract p from the total at the end.
        int32_t chk_carry;
        asm volatile(
            "sbb %%eax, %%eax\n\t"  //Essentially: Set %eax to zero then subtract carry
            "mov %%eax, %0"         //Save %eax
            : "=r" (chk_carry)      //Store 32-bit result in chk_carry (likely edx)
            : "r" (sum)             //Must happen after sum is computed
            : "%eax"
        );

        //Divide the sum by 2^64 (get most significant 64-bit word)
        result = sum >> 64;

        //Subtract p if the result is greater than p or the 128-bit sum overflowed
        if(result >= p || chk_carry == -1)
            result -= p;

    #endif

    return result;
}

//Finds the smallest positive primitive root of prime modulo using list of trial exponents.
uint64_t get_min_pr(const uint64_t p, const uint64_t* const exp_start, const uint64_t* const exp_end) {
    const uint64_t* ptr;
    uint64_t n, acc, res;

    uint64_t pi = mod_inv64(-p);
    uint64_t r2 = 1 + (uint64_t)((uint128_t) -1 % p);
    uint64_t rootp;

    //Loop through root candidates, starting at 2
    uint64_t root = 1;
    while(1) {
        //Incrementing at the start of the loop reduces the number of jump instructions
        ++root;

        //Convert root in preparation for Montgomery modular multiplication
        rootp = mm_reduce((uint128_t) root * r2, p, pi);

        //Loop through trial exponents. Back to front is more efficient as it breaks sooner
        ptr = exp_start;
        while(1) {
            n = *ptr;

            // Perform modular exponentiation by squaring, with
            //  128-bit intermediate values (faster than GMP)
            res = 1;
            acc = rootp;
            while(1) {
                if(n & 1)
                    res = mm_reduce((uint128_t) acc * res, p, pi);

                n >>= 1;
                if (n == 0) break;

                acc = mm_reduce((uint128_t) acc * acc, p, pi);
            }

            //Move on to next root candidate if result is ever 1
            if(res == 1) break;

            //Return the first root value that reaches the final exponent without breaking
            ++ptr;
            if(ptr == exp_end) return root;
        }
    }
}

//Calculates ((r ^ (p - 1) % p^2) - 1) / p, with 64-bit inputs
//Need to use GMP as there are 256-bit intermediate values
uint64_t get_residue(uint64_t root, uint64_t prime) {
    //Set g_a = root, g_b = prime^2
    mpz_set_ui(g_a, root);
    mpz_set_ui(g_b, prime);
    mpz_mul_ui(g_b, g_b, prime);

    //Calculate ((root ^ (prime - 1) % prime^2) - 1 ) / prime
    mpz_powm_ui(g_a, g_a, prime - 1, g_b);
    mpz_sub_ui(g_a, g_a, 1);               // Doing these in GMP is faster
    mpz_divexact_ui(g_a, g_a, prime);      // than exporting to a uint128_t

    //Result is necessarily 64-bit; cast down and return
    return mpz_get_ui(g_a);
}

//Thread entry point. Sieves primes, allcates exp vectors, calls the functions above
archive* find_mprs(const uint64_t chunk_start, uint32_t* count_pointer) {
    uint64_t j, modulo, root, residue, entry_index, chunk_j, vsize, fprod;
    uint64_t *exp_start, *exp_end;

    uint64_t** vectors = (uint64_t**) malloc(chunk_size * sizeof(uint64_t*));
    for(j = 0; j < chunk_size; ++j)
        vectors[j] = (uint64_t*) 1;

    //Sieve out numbers where 2n+1 is composite
    //fprintf(stderr, "Sieving primes\n"); fflush(stderr);
    prime_sieve(chunk_start, vectors);

    //Count number of entries in chunk
    uint32_t count = 0;
    for(j = 0; j < chunk_size; ++j)
        if(vectors[j])
            count++;

    //Allocate memory for list of entries before creating exponent vectors
    archive* entries = (archive*) malloc(count * sizeof(archive));

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

    //Fill vectors with (N-1)/p values for each prime number N
    //fprintf(stderr, "Factoring\n"); fflush(stderr);
    gen_exps(chunk_start, vectors);

    //Initialize 256-bit integer structs
    mpz_init(g_a);
    mpz_init(g_b);

    //Get first primitive root of each prime N
    //fprintf(stderr, "Getting MPRs\n"); fflush(stderr);
    entry_index = 0;
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

        //Calculate residue ((root^(N - 1) % N^2) - 1) / N, and print everything
        residue = get_residue(root, modulo);

        if(residue == 0)
            fprintf(stderr, "Non-generous prime found! Value = %lu (%lu)\n", modulo, root);

        entries[entry_index].prime = modulo;
        entries[entry_index].root = root;
        entries[entry_index].residue = residue;
        entry_index++;

        free(vectors[j]);
    }

    //Garbage collect 256-bit integers and vectors
    mpz_clear(g_a);
    mpz_clear(g_b);
    free(vectors);

    *count_pointer = count;
    return entries;
}

//Main entry point. Reads in start, block_size and chunk_size from argv, then splits
//the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv) {
    char* endstr_ptr;
    uint64_t sys_memory, i, j, num_chunks, max_root, maxr_prime;

    char fcsv[] = "%lu,%lu,%lu\n";
    char ftab[] = "Prime: %-19lu Root: %-4lu Residue: %lu\n";
    char* output = ftab;

    //Parse arguments, if any
    if(argc < 4) {
        puts("64-Bit Non-Generous Prime Searcher by Dylan G.\n");
        puts("Usage: ngp [c/b] start_value total_size chunk_size [> destination_file]");
        puts("    c: Format output as csv, base 10");
        puts("    b: Format output as binary data");
        puts("    (default): tabular output, base 10");
        return 0;
        
    } else if(argc > 4) {
        if(argv[1][0] == 'c')
            output = fcsv;
        else if(argv[1][0] == 'b')
            output = 0;
        else {
            fprintf(stderr, "Unrecognized option: %c\n", argv[1][0]);
            return 1;
        }
    }

    //Read in numerical inputs
    uint64_t start = strtoull(argv[argc - 3], &endstr_ptr, 10) / 2;
    uint64_t block_size = strtoull(argv[argc - 2], &endstr_ptr, 10) / 2;
    chunk_size = strtoull(argv[argc - 1], &endstr_ptr, 10) / 2;

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

    //Estimate expected number of primes in block and in (first) chunk using PNT
    uint64_t primes_in_block = noln(N_max) - noln(N_min);
    uint64_t primes_in_chunk = noln(N_min + 2 * chunk_size) - noln(N_min);

    //Estimate RAM usage. Shared: primes file & archive array. Thread: allocated vectors & ALL vector pointers.
    uint64_t shared_mem = 814000000 + primes_in_block * sizeof(archive);
    uint64_t thread_mem = primes_in_chunk * INIT_ALLOC * sizeof(uint64_t) + chunk_size * sizeof(size_t);
    uint64_t req_memory = shared_mem + thread_mem * omp_get_num_threads();

    //Check RAM capacity against prediction (assuming num_chunks isn't huge!)
    fprintf(stderr, "Predicted RAM usage: %lu MiB (capacity: %lu MiB)\n", req_memory >> 20, sys_memory >> 20);
    if(req_memory > sys_memory) {
        fprintf(stderr, "Error: insufficient system memory\n");
        return 1;
    }

    //Load in primes up through sqrt(N_max) from primes.dat
    load_primes(N_max);
    
    //Calculate the total number of chunks to be split among threads and allocate output arrays
    num_chunks = block_size / chunk_size;
    archive** lists = (archive**) malloc(num_chunks * sizeof(archive*));
    uint32_t* counts = (uint32_t*) malloc(num_chunks * sizeof(uint32_t));

    //Split chunks among threads using guided schedule since higher chunks contain fewer primes
    #pragma omp parallel for schedule(guided)
    for(i = 0; i < num_chunks; ++i)
        lists[i] = find_mprs(start + i * chunk_size, &counts[i]);
    fprintf(stderr, "Search complete, writing to output\n");

    //Format output to stdout while searching for the maximum root
    max_root = 0;
    for(i = 0; i < num_chunks; ++i) {
        if(!output)
            fwrite(lists[i], sizeof(archive), counts[i], stdout);
        else {
            for(j = 0; j < counts[i]; ++j) {
                printf(output, lists[i][j].prime, lists[i][j].root, lists[i][j].residue);
                if (lists[i][j].root < max_root) continue;

                max_root = lists[i][j].root;
                maxr_prime = lists[i][j].prime;
            }

            if(output == ftab)
                printf("\n");
        }
        free(lists[i]);
    }

    //Garbage collect
    free(lists);
    free(counts);
    free(primes);

    if(max_root)
        fprintf(stderr, "Maximum Root: %lu (%lu)\n", max_root, maxr_prime);

    return 0;
}
