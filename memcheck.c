#include <sys/sysctl.h>
#include <stdio.h>

int main() {
    int mib[2];
    size_t length;
    int64_t ram_size;

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE; // Or HW_PHYSMEM on some BSDs
    length = sizeof(ram_size);

    if (sysctl(mib, 2, &ram_size, &length, NULL, 0) == 0) {
        printf("Total RAM (bytes): %lld\n", ram_size);
        printf("Total RAM (MiB): %lld\n", ram_size / (1024 * 1024));
    } else {
        perror("Failed to get memory info via sysctl");
    }
    return 0;
}
