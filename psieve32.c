#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>

const char* output_name = "primes.dat";
unsigned int static_primes[] = {3,  5,  7, 11, 13, 17, 19, 23,
                           29, 31, 37, 41, 43, 47, 53, 59, 61,
                           67, 71, 73, 79, 83, 89, 97,101,103,
                          107,109,113,127,131,137,139,149,151,
                          157,163,167,173,179,181,191,193,197,
                          199,211,223,227,229,233,239,241,251};

int rule_out(const unsigned int* pfacs, unsigned int n, char* skip, unsigned int max)
{
    unsigned int i, j, p, count;

    for(i = 0; i < n; i++)
    {
        p = pfacs[i];

        //#pragma omp parallel for
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
    unsigned int i, j, p;

    FILE* outfile = fopen(output_name, "wb");
    if(outfile == NULL)
    {
        printf("Could not create output file %s\n", output_name);
        return 1;
    }

    char* skip = (char*) calloc(0x80000000, sizeof(char));
    skip[0] = 1;

    unsigned int total_16 = rule_out(static_primes, sizeof(static_primes)/sizeof(int), skip, 0x8000);
    unsigned int* dynamic_primes = (unsigned int*) malloc(total_16 * sizeof(int));

    i = 0;
    for(j = 0; j < 0x8000; j++)
    {
        if(skip[j]) continue;

        dynamic_primes[i] = 2 * j + 1;
        i++;
    }

    puts("Finished sieve 1, starting sieve 2"); fflush(stdout);
    unsigned int total_32 = rule_out(dynamic_primes, total_16, skip, 0x80000000);
    free(dynamic_primes);

    puts("Finished sieve 2, writing to file"); fflush(stdout);

    total_32++;
    fwrite(&total_32, sizeof(int), 1, outfile);
    p = 2;
    fwrite(&p, sizeof(int), 1, outfile);

    for(j = 0; j < 0x80000000; j++)
    {
        if(skip[j]) continue;

        p = 2 * j + 1;
        fwrite(&p, sizeof(int), 1, outfile);
    }

    free(skip);
    fclose(outfile);
    return 0;
}
