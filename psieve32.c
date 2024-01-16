#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

const char* output_name = "primes.dat";
unsigned int static_primes[54] =
    {  2,  3,  5,  7, 11, 13, 17, 19, 23,
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

        for(j = 2 * p; j < max && j >= 2 * p; j += p)
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
    unsigned int i, j;

    FILE* outfile = fopen(output_name, "wb");
    if(outfile == NULL)
    {
        printf("Could not create output file %s\n", output_name);
        return 1;
    }

    char* skip = (char*) calloc(UINT_MAX, sizeof(char));
    
    skip[0] = 1;
    skip[1] = 1;

    unsigned int total_1 = rule_out(static_primes, sizeof(static_primes)/sizeof(int), skip, USHRT_MAX);
    unsigned int* dynamic_primes = (unsigned int*) malloc(total_1 * sizeof(unsigned int));

    i = 0;
    for(j = 0; j < USHRT_MAX; j++)
    {
        if(skip[j]) continue;

        dynamic_primes[i] = j;
        i++;
    }

    puts("Finished sieve 1, starting sieve 2"); fflush(stdout);
    unsigned int total_2 = rule_out(dynamic_primes, total_1, skip, UINT_MAX);
    free(dynamic_primes);
    
    puts("Finished sieve 2, writing to file"); fflush(stdout);
    fwrite(&total_2, sizeof(int), 1, outfile);
    for(j = 0; j < UINT_MAX; j++)
    {
        if(!skip[j]){
            fwrite(&j, sizeof(int), 1, outfile);
        }
    }
        
    free(skip);
    fclose(outfile);
    return 0;
}
