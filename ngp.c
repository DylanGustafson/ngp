//Non-generous prime searching program by Dylan G.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/sysinfo.h>
#include <gmp.h>
#include <omp.h>

#define INIT_ALLOC 8 //Initial allocation: 6 prime factors (1 stdev above mean)
#define MID_ALLOC 12 //Mid allocation: 10 prime factors (3 stdev's above mean)
#define MAX_ALLOC 17 //Maximum possible (15) distinct prime factors for 64-bit value (ex: 693386350578511591)

typedef unsigned __int128 uint128_t;

//
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

/*
void print_vec(uint64_t* arr_p) {
    printf("[");
    for(size_t i = 0; i < *arr_p; i++) {
        printf("%lu", *(arr_p + i));
        if (i + 1 < *arr_p) printf(", ");
    }
    printf("]\n");
}
*/

//Returns the positive offset of the first number after 'a' divisible by 'b' (returns 0 if b divides a)
uint64_t offset(const uint64_t a, const uint64_t b)
{
    return (b - a % b) % b;
}

//Performs a sieve of eratosthenes on the chunk. Returns the number of primes needed to do so.
void prime_sieve(const uint64_t chunk_start, uint64_t** vectors)
{
    uint32_t i, p;
    uint64_t j, N_min, N_max, mstart;

    //First 2n+1 value for calculating offsets
    N_min = 2 * chunk_start + 1;
    N_max = N_min + 2 * (chunk_size - 1);

    //Loop through each prime number except 2 (since 2n+1 is always odd)
    for(i = 1; i < max_pfacs; i++)
    {
        p = primes[i];

        //Calculate offset of first N entry divisible by p.
        mstart = offset((N_min - p) >> 1, p);

        //Skip p itself if p is in the chunk
        if(N_min <= p)
            mstart += p;

        //Set false in every position divisible by p
        for(j = mstart; j < chunk_size; j += p)
            if(vectors[j])
                vectors[j] = 0;
    }
}

void gen_exps(const uint64_t chunk_start, uint64_t** vectors)
{
    uint32_t i, p;
    uint64_t j, first_pos, chunk_max, half_exp_val, newcap;
    uint128_t divisor;

    chunk_max = 2 * (chunk_start + chunk_size);

    //printf("chunk_start: %llu\tchunk_size: %llu\tchunk_max: %llu", chunk_start, chunk_size, chunk_max);

    //Primary loop: Loop through prime factors
    for(i = 0; i < max_pfacs; i++)
    {
        p = primes[i];

        //p=2 is already accounted for since the math is trivial
        if(i)
        {
            //Get the offset of the first value divisible by p.
            first_pos = offset(chunk_start, p);

            //First (2n/p)/2 value, minus 1 to fit with for loop
            half_exp_val = ((chunk_start + first_pos) / p) - 1;

            //Loop through chunk and update exp_vecs vectors using prime factor p
            for(j = first_pos; j < chunk_size; j += p)
            {
                half_exp_val++;
                if(!vectors[j]) continue;

                //update fprod (first element of each vector)
                vectors[j][2] *= p;

                //Push back prime factor value onto vector
                if (vectors[j][0] == vectors[j][1]) {
                    //printf("Reallocating for %llu:\t%llu -> %llu\n", 2*(chunk_start + j) + 1, vectors[j][1], vectors[j][1]*MULT);
                    //vectors[j] = (uint64_t*) realloc(vectors[j], vectors[j][1] * MULT * sizeof(uint64_t));
                    //vectors[j][1] *= MULT;
                    newcap = (vectors[j][1] == INIT_ALLOC) ? MID_ALLOC : MAX_ALLOC;

                    vectors[j] = (uint64_t*) realloc(vectors[j], newcap * sizeof(uint64_t));
                    vectors[j][1] = newcap;
                }

                vectors[j][vectors[j][0]] = 2 * half_exp_val;
                vectors[j][0] += 1;
            }
        }

        //Secondary Loop: Loop through powers of each prime p^e ('divisor').
        //When p=2 we have to use p^(e-1), hence the conditional start value
        for(divisor = (uint64_t) p * (i ? p : 1); divisor < chunk_max; divisor *= p)
        {
            //Get the offset of the first value divisible by p^e.
            //Should be chunk_start * 2 if chunk_start is odd
            first_pos = offset(chunk_start, divisor);

            //Break if the first number divisble by p^e is outside of the current chunk.
            if(first_pos >= chunk_size)
                break;

            //Tertiary loop: Update fprods (first element of each vector)
            for(j = first_pos; j < chunk_size; j += divisor)
                if(vectors[j])
                    vectors[j][2] *= p;
        }
    }
}

//Find the smallest positive primitive root of modulo (prime) using list of trial exponents.
uint64_t get_min_pr(const uint64_t modulo, const uint64_t* const exp_start, const uint64_t* const exp_end) {
    const uint64_t* exp_ptr;
    uint64_t n, acc, res;

    //Loop through root candidates, starting at 2
    uint64_t root = 1;
    while(1)
    {
        root++; //Incrementing at the start of the loop reduces the number of jump instructions

        //Loop through trial exponents. Back to front is more efficient as it breaks sooner
        exp_ptr = exp_start;
        while(1)
        {
            n = *exp_ptr;

            // Perform modular exponentiation by squaring, with
            //  128-bit intermediate values (faster than GMP)
            res = 1;
            acc = root;
            while(1)
            {
                // It's more efficient to just cast to a uint128 here, since %rax and
                // %rdx will already contain the 128-bit product. These casts merely
                // ensure that the modulo step uses the 128-bit subroutine __umodti3.
                if (n & 1)
                    res = ((uint128_t) acc * res) % modulo;

                n >>= 1;
                if (n == 0) break;

                acc = ((uint128_t) acc * acc) % modulo;
            }

            //Move on to next root candidate if result is ever 1
            if(res == 1) break;

            //Return the first root value that completes the secondary loop without breaking
            exp_ptr++;
            if(exp_ptr == exp_end) return root;
        }
    }
}

//Calculates ((r ^ (p - 1) % p^2) - 1) / p, with 64-bit inputs
//Need to use GMP as there are 256-bit intermediate values
uint64_t get_residue(uint64_t root, uint64_t prime)
{
    //Set g_a = root, g_b = prime^2
    mpz_set_ui(g_a, root);
    mpz_set_ui(g_b, prime);
    mpz_mul_ui(g_b, g_b, prime);

    //Calculate ((root ^ (prime - 1) % prime^2) - 1 ) / prime
    mpz_powm_ui(g_a, g_a, prime - 1, g_b);
    mpz_sub_ui(g_a, g_a, 1);               // Doing these in GMP is faster
    mpz_divexact_ui(g_a, g_a, prime);      // than exporting to a uint128

    //Result is necessarily 64-bit; cast down and return
    return mpz_get_ui(g_a);
}

//Thread entry point.
archive* find_mprs(const uint64_t chunk_start, uint32_t* count_pointer)
{
    uint32_t i, vsize;
    uint64_t j, modulo, root, residue, entry_index, fprod;
    uint64_t *exp_start, *exp_end;

    uint64_t** vectors = (uint64_t**) malloc(chunk_size * sizeof(uint64_t*));
    for(j = 0; j < chunk_size; j++)
        vectors[j] = (uint64_t*) 1;

    //Sieve out numbers where 2n+1 is composite
    //puts("Sieving primes"); fflush(stdout);
    prime_sieve(chunk_start, vectors);

    //Count number of entries in chunk
    uint32_t count = 0;
    for(j = 0; j < chunk_size; j++)
        if(vectors[j])
            count++;

    //Allocate memory for list of entries before creating exponent vectors
    archive* entries = (archive*) malloc(count * sizeof(archive));

    //Create array of vectors which contain test exponent values for each prime
    //std::vector<uint64_t>* exp_vecs = new std::vector<uint64_t>[chunk_size];

    for(j = 0; j < chunk_size; j++) {
        if(!vectors[j])
            continue;

        vectors[j] = (uint64_t*) malloc(INIT_ALLOC * sizeof(int64_t));
        vectors[j][0] = 3;
        vectors[j][1] = INIT_ALLOC;
        vectors[j][2] = 1;

        //printf("Allocated vector %llu: %llu --> %llu, %llu\n", (uint64_t)i, (uint64_t)arr_pp[i], arr_pp[i][0], arr_pp[i][0]);
    }

    //Fill vectors with (N-1)/p values for each prime number N
    //puts("Factoring"); fflush(stdout);
    gen_exps(chunk_start, vectors);

    //Initialize 256-bit integer structs
    mpz_init(g_a);
    mpz_init(g_b);

    //puts("Getting MPRs"); fflush(stdout);
    //Get first primitive root of each prime N
    entry_index = 0;
    for(j = 0; j < chunk_size; j++)
    {
        if(!vectors[j]) continue;

        //print_vec(vectors[j]);

        //Calculate the test prime N
        modulo = 2 * (chunk_start + j) + 1;
        vsize = vectors[j][0];

        //fprod is just N-1, meaning all prime factors of N-1 are less than sqrt(N-1)
        //In this case, simply overwrite fprod with the first exponent, (N-1)/2

        exp_start = vectors[j] + 2;
        exp_end = vectors[j] + vsize;

        if(vectors[j][2] != chunk_start + j)
        {
            fprod = vectors[j][2] * 2;

            if(vectors[j][1] > vsize) {
                vectors[j][vsize] = fprod;
                vectors[j][2] = chunk_start + j;
                exp_end += 1;
            }
            else
            {
                for(i = 2; i < vsize - 2; i++)
                    vectors[j][i] = vectors[j][i + 1];

                vectors[j][vsize - 1] = fprod;
                vectors[j][1] = chunk_start + j;
                exp_start -= 1;
            }
        }

        //Find first primitive root of N
        root = get_min_pr(modulo, exp_start, exp_end);

        //Calculate residue ((root^(N - 1) % N^2) - 1) / N, and print everything
        residue = get_residue(root, modulo);

        if(residue == 0)
            fprintf(stderr, "Non-generous prime found! Value = %lu\n", modulo);

        entries[entry_index].prime = modulo;
        entries[entry_index].root = root;
        entries[entry_index].residue = residue;
        entry_index++;

        free(vectors[j]);
    }

    //Garbage collect 256-bit integers
    mpz_clear(g_a);
    mpz_clear(g_b);

    //Garbage collect and return
    free(vectors);

    *count_pointer = count;
    return entries;
}

void load_primes(uint64_t Nmax)
{
    uint32_t i, pfile_count;
    size_t fread_out;

    //Open primes.dat to get primes needed for sieve
    FILE* pfile = fopen(primes_file, "rb");
    if(pfile == NULL)
    {
        fprintf(stderr, "Could not find \"%s\"\n", primes_file);
        exit(1);
    }

    //First 4 bytes is total number of primes
    if(fread(&pfile_count, sizeof(uint32_t), 1, pfile) != 1) {
        fprintf(stderr, "Format error in \"%s\"\n", primes_file);
        exit(1);
    }
    
    //The rest of the file is loaded straight into the array in 4KB increments
    primes = (uint32_t*) malloc(pfile_count * sizeof(uint32_t));
    max_pfacs = 0;
    while(max_pfacs < pfile_count) {
        max_pfacs += fread(&primes[max_pfacs], sizeof(uint32_t), 1024, pfile);
        
        //Break when p^2 is larger than the largest value in the block
        if((uint64_t) primes[max_pfacs - 1] * primes[max_pfacs - 1] > Nmax)
            break;
    }
    fclose(pfile);

    //Truncate primes array (if possible) to free up memory for next steps
    primes = (uint32_t*) realloc(primes, max_pfacs * sizeof(uint32_t));
}

//Main entry point. Reads in upper limit and chunk_size from argv, then splits
//the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv)
{
    //puts("Initializing");
    char* endstr_ptr;
    uint64_t j;
    uint32_t num_chunks;

    char fcsv[] = "%lu,%lu,%lu\n";
    char fstdout[] = "Prime: %19lu    Root:%4lu    Residue: %lu\n";
    char* output = fstdout;

    if(argc < 4) {
        puts("64-Bit Non-Generous Prime Searcher by Dylan G.\n");
        puts("Usage: ngp [-cbs] start_value total_size chunk_size");
        puts("   -c: Format output as csv (base 10)");
        puts("// -b: Format output as binary");
        puts("// -s: Silent search, no output to stdout");
        return 0;
        
    } else if(argc > 4 && argv[1][1] == 'c') {
        output = fcsv;
    }

    //Read in bounds and chunk size from command line (defaults to 0)
    uint64_t start = strtoull(argv[argc - 3], &endstr_ptr, 10) / 2;
    uint64_t block_size = strtoull(argv[argc - 2], &endstr_ptr, 10) / 2;

    //Load in primes up through sqrt(block end) from primes.dat
    load_primes(2 * (start + block_size));

    chunk_size = strtoull(argv[argc - 1], &endstr_ptr, 10) / 2;
    num_chunks = block_size / chunk_size;
    archive** lists = (archive**) malloc(num_chunks * sizeof(archive*));
    uint32_t* counts = (uint32_t*) malloc(num_chunks * sizeof(uint32_t));
    
    struct sysinfo info;
    if (sysinfo(&info) != 0) {
        perror("sysinfo");
        return 1;
    }
    
    //Check RAM size against rough prediction of usage
    if(chunk_size * INIT_ALLOC + 1e9 > info.totalram / omp_get_num_threads()) {
        fprintf(stderr, "Error - insufficient system memory. Try a smaller chunk size.\n");
        return 1;
    }

    //Use OpenMP to split chunks among threads
    #pragma omp parallel for
    for(uint32_t i = 0; i < num_chunks; i++)
        lists[i] = find_mprs(start + i * chunk_size, &counts[i]);

    //puts("Displaying"); fflush(stdout);
    for(uint32_t i = 0; i < num_chunks; i++)
    {
        for(j = 0; j < counts[i]; j++)
            printf(output, lists[i][j].prime, lists[i][j].root, lists[i][j].residue);
        printf("\n");
        free(lists[i]);
    }
    free(lists);
    free(counts);

    //Dont need list of primes anymore
    free(primes);
    return 0;
}
