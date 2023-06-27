#ifndef FMD_PERMUTATION_PARSER_H
#define FMD_PERMUTATION_PARSER_H
#include "fmd_vector.h"
#include <stdio.h>

static fmd_vector_t *permutation_parse(FILE* file) {
    fmd_vector_t *vec;
    fmd_vector_init(&vec, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        int_t val;
        val = strtoll(line, NULL, 10);
        fmd_vector_append(vec, (void*)val);
    }
    fmd_vector_fit(vec);
    return vec;
}

#endif //FMD_PERMUTATION_PARSER_H
