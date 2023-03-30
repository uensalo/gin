#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/fmd_bitstream.h"

struct timespec t1,t2;

/*
 *
 *
clock_gettime(CLOCK_REALTIME, &t1);
clock_gettime(CLOCK_REALTIME, &t2);
double time1 = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec)*(1e-9);
 */

int main() {
    fmd_bs_t *bs;
    int no_writes = 2022;
    fmd_bs_init_reserve(&bs,1);


    int pointer = 0;
    for (int i = 0; i < no_writes; i++) {
        word_t val = (i % 7) * (i % 11) * (i % 17);
        word_t sz = (i * i) % 13;
        fmd_bs_init_write_word(bs, pointer, val, sz);
        printf("%d ", ~mask_shift_left_64[sz] & val);
        pointer += sz;
    }
    printf("\n");
    pointer = 0;
    for (int i = 0; i < no_writes; i++) {
        word_t val = 0;
        word_t sz = (i * i) % 13;
        fmd_bs_init_read_word(bs, pointer, sz, &val);
        printf("%d ", val);
        pointer += sz;
    }

}
