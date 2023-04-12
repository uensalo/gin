#include <stdio.h>
#include <stdlib.h>
#include <fmd_table.h>
#include <fmd_vector.h>
#include <fmd_bitstream.h>
#include <time.h>
#include <fmd_graph.h>
#include <fmd_fmi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

/*
// fmi unit test
fmd_string_t *str,*query;
fmd_string_init_cstr(&str, "abracadabra");
fmd_string_init_cstr(&query, "q");

fmd_fmi_t *fmi;
fmd_fmi_init(&fmi, str, 5, 5);

for(int i = 0; i < fmi->alphabet->size; i++) {
    void* encoding;
    char c = (char)fmi->alphabet->data[i];
    fmd_table_lookup(fmi->e2c, c, &encoding);
    printf("%c %d %d\n", fmi->alphabet->data[i], encoding, fmi->char_counts[i]);
}

for(int j = 0; j < fmi->no_chars; j++) {
    int encoding = fmd_fmi_get(fmi, j);
    printf("%d ", encoding);
}
printf("\n");


for(int j = 0; j < fmi->no_chars; j++) {
    for (int i = 1; i < fmi->alphabet->size; i++) {
        void *encoding;
        char c = (char) fmi->alphabet->data[i];
        fmd_table_lookup(fmi->c2e, c, &encoding);
        pos_t rank = fmd_fmi_rank(fmi, encoding, j);
        printf("%d ", rank);
    }
    printf("\n");
}

count_t count = fmd_fmi_query_count(fmi, query);
printf("Count: %llu\n", count);

fmd_vector_t *vec = fmd_fmi_query_locate(fmi, query);
fmd_vector_sort(vec);
printf("Locations: ");
for(int i = 0; i < vec->size; i++) {
    printf("%d ", vec->data[i]);
}
printf("\n");
 */

#include "fmd_interval_merge_tree.h"

void print_intervals(fmd_vector_t *intervals) {
    for (int_t i = 0; i < intervals->size; i++) {
        fmd_imt_interval_t *interval = intervals->data[i];
        printf("(%ld, %ld)", interval->lo, interval->hi);
        if (i < intervals->size - 1) {
            printf(", ");
        }
    }
    printf("\n");
}

int main() {
    fmd_vector_t *intervals1, *intervals2;
    fmd_vector_init(&intervals1, 10, &fmd_fstruct_imt_interval);
    fmd_vector_init(&intervals2, 10, &fmd_fstruct_imt_interval);

    fmd_imt_interval_t i1_1 = {1, 3};
    fmd_imt_interval_t i1_2 = {7, 9};
    fmd_imt_interval_t i1_3 = {15, 18};

    fmd_imt_interval_t i2_1 = {5, 6};
    fmd_imt_interval_t i2_2 = {11, 14};
    fmd_imt_interval_t i2_3 = {19, 20};

    fmd_vector_append(intervals1, &i1_1);
    fmd_vector_append(intervals1, &i1_2);
    fmd_vector_append(intervals1, &i1_3);

    fmd_vector_append(intervals2, &i2_1);
    fmd_vector_append(intervals2, &i2_2);
    fmd_vector_append(intervals2, &i2_3);

    printf("Intervals 1: ");
    print_intervals(intervals1);
    printf("Intervals 2: ");
    print_intervals(intervals2);

    fmd_vector_t *merged_intervals = fmd_imt_merge_intervals(intervals1, intervals2);
    printf("Merged Intervals: ");
    print_intervals(merged_intervals);

    fmd_vector_free(intervals1);
    fmd_vector_free(intervals2);
    fmd_vector_free(merged_intervals);

    return 0;
}