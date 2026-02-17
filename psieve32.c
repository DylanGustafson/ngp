//Dylan G.
//Generates all primes below 2^32 and store them in binary file
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

//32768 * 15. Splits 0x8000 to 0x80000000 sieve into 4369 chunks
#define CHUNK_SIZE 491520

typedef uint32_t u32;

// All primes below 2^8
const u32 primes8[] = {
      2,  3,  5,  7, 11, 13, 17, 19, 23,
     29, 31, 37, 41, 43, 47, 53, 59, 61,
     67, 71, 73, 79, 83, 89, 97,101,103,
    107,109,113,127,131,137,139,149,151,
    157,163,167,173,179,181,191,193,197,
    199,211,223,227,229,233,239,241,251
};

// Pattern for 3, 5, and 7
const char pattern[] = {
    0,1,1,0,1,0,0,0,1,0,0,1,0,1,0,
    0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,
    0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
    0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,
    0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,
    0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,
    0,0,1,0,1,0,0,1,0,0,0,1,0,1,1
};

const u32 p8count = sizeof(primes8)/sizeof(u32);
const u32 pattern_size = sizeof(pattern);
const char* output_name = "primes.dat";
char* keep;

// Performs a Sieve of Eratosthenes on a chunk of keep[] array using pfacs
// Returns the total prime count in chunk
uint32_t rule_out(const u32* pfacs, const u32 npfacs, const u32 chunk_start, const u32 chunk_stop)
{
    u32 i, j, p, mstart, count;

    // Copy in the tail end of the pattern up to the first multiple of 105
    u32 copy_start = (chunk_start + (pattern_size + 1) / 2) % pattern_size;
    if (__builtin_expect(copy_start != 0, 1)) {
        memcpy(&keep[chunk_start], &pattern[copy_start], pattern_size - copy_start);
        copy_start = pattern_size - copy_start;
    }

    // Copy in the full pattern until it can't fit
    for (j = chunk_start + copy_start; j + pattern_size <= chunk_stop; j += pattern_size)
        memcpy(&keep[j], pattern, pattern_size);

    // Copy in the beginning of the pattern up into the last bit of the chunk
    u32 tail_size = chunk_stop - j;
    if (__builtin_expect(tail_size != 0, 1))
        memcpy(&keep[j], pattern, tail_size);

    // Loop through given primes for sieve
    for(i = 4; i < npfacs; i++) {
        p = pfacs[i];
        mstart = p * p / 2;

        // Otherwise just get the offset of the first N divisible by p
        if (mstart < chunk_start)
            mstart = chunk_start + p - 1 - (chunk_start - mstart - 1) % p;

        // Avoid unnecessary compiler optimization branch
        if (p == 1) __builtin_unreachable();

        // Zero every position divisible by p
        for (j = mstart; j < chunk_stop; j += p)
            if (keep[j])
                keep[j] = 0;
    }

    // Count up primes in list
    count = 0;
    for (j = chunk_start; j < chunk_stop; ++j)
        if (keep[j])
            ++count;

    return count;
}

int main(int argc, char* argv[])
{
    u32 count, j, p;

    FILE* outfile = fopen(output_name, "wb");
    if(outfile == NULL) {
        printf("Could not create output file %s\n", output_name);
        return 1;
    }

    keep = (char*) aligned_alloc(64, 0x80000000);

    // Use primes below 2^8 to get all primes below 2^16
    u32 p16count = p8count + rule_out(primes8, p8count, 0x80, 0x8000);

    // Create and fill array of primes below 2^16
    u32* primes16 = (u32*) malloc(p16count * sizeof(u32));
    memcpy(primes16, primes8, sizeof(primes8));

    count = p8count;
    for (j = 0x80; j < 0x8000; ++j) {
        if(!keep[j]) continue;

        primes16[count] = 2 * j + 1;
        ++count;
    }

    // Use primes below 2^16 to get all primes below 2^32
    puts("Starting 32-bit sieve"); fflush(stdout);

    u32 num_chunks = (0x80000000 - 0x8000) / CHUNK_SIZE;
    u32* chunk_totals = (u32*) malloc(num_chunks * sizeof(u32));
    u32** primes32 = (u32**) malloc(num_chunks * sizeof(u32*));

    #pragma omp parallel for
    for (u32 i = 0; i < num_chunks; ++i) {
        u32 start = 0x8000 + CHUNK_SIZE * i;
        u32 stop = start + CHUNK_SIZE;

        chunk_totals[i] = rule_out(primes16, p16count, start, stop);
        primes32[i] = (u32*) malloc(chunk_totals[i] * sizeof(u32));

        u32 p32count = 0;
        for (u32 k = start; k < stop; ++k) {
            if (!keep[k]) continue;

            primes32[i][p32count] = 2 * k + 1;
            ++p32count;
        }
    }

    puts("Finished 32-bit sieve, writing to file"); fflush(stdout);

    // Get grand total for 32-bit primes
    u32 p32count = p16count;
    for (j = 0; j < num_chunks; ++j)
        p32count += chunk_totals[j];

    // Write total number of primes to first dword
    fwrite(&p32count, sizeof(u32), 1, outfile);

    // Write
    fwrite(primes16, sizeof(u32), p16count, outfile);
    for (j = 0; j < num_chunks; ++j) {
        fwrite(primes32[j], sizeof(u32), chunk_totals[j], outfile);
        free(primes32[j]);
    }

    free(primes32);
    free(chunk_totals);
    free(primes16);
    free(keep);
    fclose(outfile);
    return 0;
}
