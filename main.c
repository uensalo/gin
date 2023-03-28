#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/fmd_bitstream.h"

struct timespec t1,t2;

int main() {
    fmd_bs_t *bs;
    int no_writes = 100;
    int seed = 4325534;
    fmd_bs_init_reserve(&bs,10000);

    srand(seed);
    clock_gettime(CLOCK_REALTIME, &t1);
    int pointer = 0;
    for (int i = 0; i < no_writes; i++) {
        //printf("%d\n", i);
        word_t val = (word_t)rand();
        word_t sz = rand() % 32;
        fmd_bs_init_write_word(bs, pointer, val, sz);
        pointer += sz;
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    double time1 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec)*(1e-9);
    pointer = 0;
    srand(seed);
    for (int i = 0; i < no_writes; i++) {
        //printf("%d\n", i);
        word_t val=0;
        word_t valr = (word_t)rand();
        word_t sz = rand() % 32;
        valr = valr & (0xFFFFFFFFFFFFFFFF >> (64-sz));
        fmd_bs_init_read_word(bs, pointer, sz, &val);
        if(valr != val) {
            printf("valr %d vs val %d at %d\n", valr, val, i);
            exit(-1);
        }
        pointer += sz;
    }


    printf("write  time: %lf\n", time1);
}
