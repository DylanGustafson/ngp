//Dylan G.
//Generates all primes below 2^32 and store them in binary file
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define BLOCK_SIZE 0x8000

const char* output_name = "primes.dat";

//All (odd) primes below 2^8
uint32_t static_primes[] = {  3,  5,  7, 11, 13, 17, 19, 23,
                         29, 31, 37, 41, 43, 47, 53, 59, 61,
                         67, 71, 73, 79, 83, 89, 97,101,103,
                        107,109,113,127,131,137,139,149,151,
                        157,163,167,173,179,181,191,193,197,
                        199,211,223,227,229,233,239,241,251};

//Perform a Sieve of Eratosthenes using primes below 2^m to get primes below 2^(2m)
uint32_t rule_out(const uint32_t* pfacs, const uint32_t n, char* skip, const uint32_t max)
{
    uint32_t i, j, p, count;

    //Loop through given primes for sieve
    for(i = 0; i < n; i++)
    {
        p = pfacs[i];

        //Start from p^2 and flag all subsequent mutiples of p
        for(j = (p * p - 1) / 2; j < max; j += p)
            skip[j] = 1;
    }

    //Count up primes in list
    count = 0;
    for(j = 0; j < max; j++)
        if(!skip[j])
            count++;

    return count;
}

int main(int argc, char* argv[])
{
    uint32_t pblock[BLOCK_SIZE];
    uint32_t count, j, p;

    FILE* outfile = fopen(output_name, "wb");
    if(outfile == NULL)
    {
        printf("Could not create output file %s\n", output_name);
        return 1;
    }

    char* skip = (char*) calloc(0x80000000, sizeof(char));
    skip[0] = 1;

    //Use primes below 2^8 to get all primes below 2^16
    uint32_t total_16 = rule_out(static_primes, sizeof(static_primes)/sizeof(uint32_t), skip, 0x8000);

    //Create and fill array of primes below 2^16
    uint32_t* dynamic_primes = (uint32_t*) malloc(total_16 * sizeof(uint32_t));
    count = 0;
    for(j = 0; j < 0x8000; j++)
    {
        if(skip[j]) continue;

        dynamic_primes[count] = 2 * j + 1;
        count++;
    }

    //Use primes below 2^16 to get all primes below 2^32
    puts("Starting 32-bit sieve"); fflush(stdout);
    uint32_t total_32 = rule_out(dynamic_primes, total_16, skip, 0x80000000);
    free(dynamic_primes);
    puts("Finished 32-bit sieve, writing to file"); fflush(stdout);

    //Increment to account for p=2, then write total count to start of file
    total_32++;
    fwrite(&total_32, sizeof(uint32_t), 1, outfile);

    //Hard code p=2. Not strictly necessary, but the output file may be used for other projects
    p = 2;
    fwrite(&p, sizeof(uint32_t), 1, outfile);

    //Fill up static pblock array with primes and flush to outfile when it hits BLOCK_SIZE
    count = 0;
    for(j = 0; j < 0x80000000; j++)
    {
        if(skip[j]) continue;

        pblock[count] = 2 * j + 1;
        count++;
        if(count < BLOCK_SIZE) continue;

        fwrite(pblock, sizeof(uint32_t), BLOCK_SIZE, outfile);
        count = 0;
    }

    //Write last pblock since the loop above ended on a composite number (2^64-1)
    fwrite(pblock, sizeof(uint32_t), count, outfile);

    free(skip);
    fclose(outfile);
    return 0;
}
