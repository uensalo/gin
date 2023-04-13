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

/*
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
 */

/*
#include <stdio.h>
#include <assert.h>

void test_fmd_vector_argsort() {
    // Initialize test data and the fmd_vector_t instance
    int_t test_data[] = {141, 53, 177, 17, 221, 5, 201, 234, 241, 239, 69, 162, 245, 246, 102, 88, 218, 116, 57, 214, 173, 210, 142, 128, 215, 99, 131, 50, 196, 233, 86, 187, 149, 160, 100, 89, 32, 72, 112, 39, 123, 22, 109, 25, 227, 40, 97, 93, 146, 192, 182, 133, 70, 144, 83, 178, 157, 169, 30, 202, 79, 20, 217, 37, 84, 62, 103, 51, 24, 175, 8, 13, 165, 28, 248, 190, 212, 61, 56, 126, 240, 198, 226, 71, 207, 135, 191, 67, 23, 189, 154, 147, 220, 166, 132, 10, 29, 199, 219, 249, 152, 104, 213, 92, 7, 193, 18, 75, 129, 31, 47, 45, 106, 42, 181, 159, 73, 164, 127, 44, 120, 232, 48, 34, 9, 243, 101, 60, 3, 204, 151, 155, 222, 15, 172, 55, 87, 12, 231, 94, 153, 38, 184, 206, 238, 4, 209, 122, 203, 161, 11, 171, 137, 54, 46, 81, 35, 95, 96, 247, 6, 14, 168, 2, 19, 76, 179, 194, 170, 236, 90, 139, 124, 134, 180, 16, 41, 108, 43, 98, 66, 107, 111, 64, 21, 143, 174, 156, 176, 26, 130, 27, 78, 251, 150, 113, 117, 49, 52, 121, 205, 225, 93, 138, 91, 237, 33, 85, 77, 110, 250, 195, 36, 148, 229, 145, 1, 208, 63, 200, 58, 244, 0, 59, 68, 228, 197, 223, 235, 167, 65, 185, 163, 186, 118, 115, 230, 74, 119, 188, 211, 80, 216, 158, 183, 242, 253};
    size_t n = sizeof(test_data) / sizeof(test_data[0]);
    fmd_vector_t *vec;
    fmd_vector_init(&vec, n, &prm_fstruct); // Replace `some_fstruct` with the appropriate fstruct

    // Fill the vector with the test data
    for (size_t i = 0; i < n; i++) {
        fmd_vector_append(vec, test_data[i]);
    }

    fmd_vector_t *original = fmd_vector_copy(vec);

    // Perform argsort
    fmd_vector_t *args_sorted;
    fmd_vector_argsort(&args_sorted, vec);

    // Verify the sorted indices are correct
    for (size_t i = 0; i < args_sorted->size; i++) {
        //printf("%d, %d\n", (int_t)vec->data[i], original->data[i]);
        assert((int_t)vec->data[i] == (int_t)original->data[(int_t)args_sorted->data[i]]);
    }

    // Clean up
    fmd_vector_free(vec);
    fmd_vector_free(args_sorted);
}

int main() {
    test_fmd_vector_argsort();
    printf("All tests passed.\n");
    return 0;
}
 */
#include <stdio.h>
#include <math.h>

#include "fmd_fmd.h"

int main() {
    fmd_fmd_init_pcodes_fixed_binary_helper('1', '0', 100);
    return 0;
}
