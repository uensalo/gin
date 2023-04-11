#include <stdio.h>
#include <stdlib.h>
#include <fmd_table.h>
#include <fmd_vector.h>
#include <fmd_bitstream.h>
#include <time.h>
#include <fmd_graph.h>
#include <fmd_fmi.h>

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

int main() {
    return 0;
}