//Dylan G.
//Generates all primes below 2^32 and store them in binary
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

const char* output_name = "primes.dat";
uint32_t static_primes[] = {  3,  5,  7, 11, 13, 17, 19, 23,
                         29, 31, 37, 41, 43, 47, 53, 59, 61,
                         67, 71, 73, 79, 83, 89, 97,101,103,
                        107,109,113,127,131,137,139,149,151,
                        157,163,167,173,179,181,191,193,197,
                        199,211,223,227,229,233,239,241,251};

uint32_t rule_out(const uint32_t* pfacs, const uint32_t n, char* skip, const uint32_t max)
{
    uint32_t i, j, p, count;

    for(i = 0; i < n; i++)
    {
        p = pfacs[i];

        for(j = (3 * p - 1) / 2; j < max; j += p)
            skip[j] = 1;
    }

    count = 0;
    for(j = 0; j < max; j++)
        if(!skip[j])
            count++;

    return count;
}

int main(int argc, char* argv[])
{
    uint32_t i, j, p;

    FILE* outfile = fopen(output_name, "wb");
    if(outfile == NULL)
    {
        printf("Could not create output file %s\n", output_name);
        return 1;
    }

    char* skip = (char*) calloc(0x80000000, sizeof(char));
    skip[0] = 1;

    uint32_t total_16 = rule_out(static_primes, sizeof(static_primes)/sizeof(uint32_t), skip, 0x8000);
    uint32_t* dynamic_primes = (uint32_t*) malloc(total_16 * sizeof(uint32_t));

    i = 0;
    for(j = 0; j < 0x8000; j++)
    {
        if(skip[j]) continue;

        dynamic_primes[i] = 2 * j + 1;
        i++;
    }

    puts("Finished sieve 1, starting sieve 2"); fflush(stdout);
    uint32_t total_32 = rule_out(dynamic_primes, total_16, skip, 0x80000000);
    free(dynamic_primes);
    
    puts("Finished sieve 2, writing to file"); fflush(stdout);
    total_32++;
    fwrite(&total_32, sizeof(uint32_t), 1, outfile);
    
    p = 2;
    fwrite(&p, sizeof(uint32_t), 1, outfile);

    for(j = 0; j < 0x80000000; j++)
    {
        if(skip[j]) continue;

        p = 2 * j + 1;
        fwrite(&p, sizeof(uint32_t), 1, outfile);
    }

    free(skip);
    fclose(outfile);
    return 0;
}
