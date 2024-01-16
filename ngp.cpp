//Queneau searching program by Dylan G.
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gmp.h>
//#include <omp.h>

//IO file names
const char* primes_file = "primes.dat";

//Global vars that remain constant across each thread
long chunk_size;
int* primes;
int max_pfacs;

//Returns the positive offset of the first number after 'a' divisible by 'b' (returns 0 if b divides a)
long offset(const long a, const long b)
{
    return (b - a % b) % b;
}

//Calculates [a^n mod b] in O(log n) time
__int128 powmod(__int128 a, long n, const long b)
{
    __int128 res = 1;
    while (n > 0)
    {
        if (n & 1)
            res = res * a % b;
        a = a * a % b;
        n >>= 1;
    }
    return res;
}

//Calculates ((a ^ (p - 1) % p^2) - 1) / p, with 64-bit inputs
//Need to use GMP as there are 256-bit intermediate values
unsigned long get_residue(unsigned long a, unsigned long p)
{
    mpz_t base, mod, r1, r2;

    //Allocate mpz structs
    mpz_init(base);
    mpz_init(mod);
    mpz_init(r1);
    mpz_init(r2);

    //Set base = a, mod = p^2
    mpz_set_ui(base, a);
    mpz_set_ui(r1, p);
    mpz_mul_ui(mod, r1, p);

    //Calculate ((base^(p-1) % mod)-1)/p and cast result to 64-bit
    mpz_powm_ui(r1, base, p - 1, mod);
    mpz_sub_ui(r2, r1, 1);
    mpz_divexact_ui(r1, r2, p);
    a = mpz_get_ui(r1);

    //garbage collect
    mpz_clear(r2);
    mpz_clear(r1);
    mpz_clear(mod);
    mpz_clear(base);

    return a;
}

//Performs a sieve of eratosthenes on the chunk. Returns the number of primes needed to do so.
void prime_sieve(const long chunk_start, char *candidates)
{
    int i, p;
    long j, N_min, N_max, mstart;

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
            if(candidates[j])
                candidates[j] = 0;
    }
}

//Get the maximum number of distinct prime factors a number below Nmax can have, minus one
int get_max_vsize(const long Nmax)
{
    long product = 1;
    for(int i = 0; i < max_pfacs; i++)
    {
        product *= primes[i];
        if(product >= Nmax) return i - 1;
    }

    printf("This error should never happen\n");
    exit(1);
}

void gen_exps(const long chunk_start, char* candidates, std::vector<long>* exp_vals)
{
    int i, p;
    long j, first_pos, step, chunk_stop, order, pm, count, half_exp_val;
    __int128 divisor;

    for(j = 0; j < chunk_size; j++)
        if(candidates[j])
            exp_vals[j].push_back(2);

    chunk_stop = 2 * (chunk_start + chunk_size);

    //Primary loop: Loop through prime factors
    for(i = max_pfacs - 1; i >= 0; i--)
    {
        p = primes[i];

        //p=2 is already accounted for since the math is trivial
        if(i)
        {
            //Get the offset of the first value divisible by p.
            first_pos = offset(chunk_start, p);

            //First (2n/p)/2 value, minus 1 to fit with for loop
            half_exp_val = ((chunk_start + first_pos) / p) - 1;

            //Loop through chunk and update exp_vals vectors using prime factor p
            for(j = first_pos; j < chunk_size; j += p)
            {
                half_exp_val++;
                if(!candidates[j]) continue;

                //update fprod (first element of each vector)
                exp_vals[j][0] *= p;

                //Add exp_value to vector
                exp_vals[j].push_back(2 * half_exp_val);
            }
        }

        //Secondary Loop: Loop through powers of each prime p^e ('divisor').
        //When p=2 we have to use p^(e-1), hence the conditional start value
        for(divisor = (long)p * (i ? p : 1); divisor < chunk_stop; divisor *= p)
        {
            //Get the offset of the first value divisible by p^e.
            //Should be chunk_start * 2 if chunk_start is odd
            first_pos = offset(chunk_start, divisor);

            //Break if the first number divisble by p^e is outside of the current chunk.
            if(first_pos >= chunk_size)
                break;

            //Tertiary loop: Update fprods (first element of each vector)
            for(j = first_pos; j < chunk_size; j += divisor)
            {
                if(candidates[j])
                    exp_vals[j][0] *= p;
            }
        }
    }
}

char chk_pr(const long root, const long modulo, const long *exp_diffs, const int vsize)
{
    __int128 res_prev = 1;
    for(int i = 0; i < vsize; i++)
    {
        res_prev = (res_prev * powmod(root, exp_diffs[i], modulo)) % modulo;
        if(res_prev == 1)
            return 0;
    }

    return 1;
}

//Thread entry point.
void find_mprs(const long chunk_start)
{
    int num_pfacs;
    int i, p, vsize, max_vsize;
    long j, exp_prev, modulo, root, a;
    int vstart;
    long* exp_diffs;
    long residue;

    long pfsum = 0;
    long counter = 0;

    max_vsize = get_max_vsize(2 * (chunk_start + chunk_size));
    exp_diffs = (long*) malloc((max_vsize + 1) * sizeof(long));

    char* candidates = (char*) malloc(chunk_size);
    for(j = 0; j < chunk_size; j++)
        candidates[j] = 1;

    //Sieve out numbers where 2n+1 is composite, and update num_pfacs
    puts("Sieving primes"); fflush(stdout);
    prime_sieve(chunk_start, candidates);

    //Create array of vectors which contain test exponent values for each prime
    std::vector<long>* exp_vals = new std::vector<long>[chunk_size];

    //Fill vectors with (N-1)/p values for each prime number N
    puts("Factoring"); fflush(stdout);
    gen_exps(chunk_start, candidates, exp_vals);

    puts("Getting MPRs"); fflush(stdout);
    //Get first primitive root of each prime N
    for(j = 0; j < chunk_size; j++)
    {
        if(!candidates[j]) continue;

        //Number of distinct prime factors (minus 1 when N-1 is unusual)
        //p=2 case is not stored, however factor prod is always stored
        vsize = exp_vals[j].size();

        pfsum += vsize;
        counter++;

        //Calculate the test prime N
        modulo = 2 * (chunk_start + j) + 1;

        //Skip factor product in first slot unless N-1 is unusual
        vstart = 1;
        if(exp_vals[j][0] != modulo - 1)
        {
            vsize++;
            vstart = 0;
        }

        //Create vector of differences between (N-1)/p value to make powmod faster
        exp_prev = 0;
        for(i = 0; i < vsize - 1; i++)
        {
            exp_diffs[i] = exp_vals[j][vstart + i] - exp_prev;
            exp_prev = exp_vals[j][vstart + i];
        }

        //Last (N-1)/p value is when p=2, which wasnt actually stored since its just n
        exp_diffs[vsize - 1] = chunk_start + j - exp_prev;

        exp_prev = 0;
        printf("\n%li: ", modulo);
        for(i = 0; i < vsize; i++)
        {
            exp_prev += exp_diffs[i];
            printf("%li ", exp_prev);
        }

        //Find first primitive root of N
        for(root = 2; root < modulo; root++)
            if(chk_pr(root, modulo, exp_diffs, vsize))
                break;

        //Calculate residue ((root^(N - 1) % N^2) - 1) / N, and print everything
        residue = get_residue(root, modulo);
        //printf("Prime: %li\t\tRoot: %3li\tA: %li\n", modulo, root, residue);fflush(stdout);
    }

    printf("\nAvg N pfacs: %f\tMax N pfacs: %li\n", (double)pfsum / counter, max_vsize);

    //Garbage collect and return
    delete[] exp_vals;
    free(candidates);
    free(exp_diffs);
}

void load_primes(long Nmax)
{
    int i, pfile_count;

    //Open primes.dat to get primes needed for sieve
    FILE* pfile = fopen(primes_file, "rb");
    if(pfile == NULL)
    {
        printf("Could not find \"%s\"\n", primes_file);
        exit(1);
    }

    //First 4 bytes is total number of primes
    fread(&pfile_count, sizeof(int), 1, pfile);
    primes = (int*) malloc(pfile_count * sizeof(int));

    for(i = 0; i < pfile_count; i++)
    {
        fread(&primes[i], sizeof(int), 1, pfile);
        if((long)primes[i] * primes[i] > Nmax)
            break;
    }

    fclose(pfile);
    if(i >= pfile_count)
    {
        printf("Error: Not enough primes in %s\n", primes_file);
        free(primes);
        exit(1);
    }

    max_pfacs = i;
    primes = (int*) realloc(primes, max_pfacs * sizeof(int));
}

//Main entry point. Reads in upper limit and chunk_size from argv, then splits
//the total range by chunk_size to be divided amongst parallel threads using OpenMP.
int main(int argc, char **argv)
{
    puts("Initializing");
    long j, total_count;

    //Read in bounds and chunk size from command line (defaults to 0)
    long start = atol(argv[argc - 3]) / 2;
    long block_size = atol(argv[argc - 2]) / 2;
    chunk_size = atol(argv[argc - 1]) / 2;


    // For now: single thread
    chunk_size = block_size;


    //Load in primes up through sqrt(block end) from primes.dat
    load_primes(2 * (start + block_size));

    //Use OpenMP to split chunks among threads
    //#pragma omp parallel for
    //for(int i = 0; i < nchunks; i++)
    find_mprs(start);

    //Dont need list of primes anymore
    free(primes);
    return 0;
}
