#ifndef FMD_FMD_VECTOR_H
#define FMD_FMD_VECTOR_H
#include "fmd_common.h"

#define FMD_VECTOR_INIT_SIZE 8

typedef struct fmd_vector_ {
    void **data;
    int_t size;
    int_t capacity;
    fmd_fstruct_t *f;
} fmd_vector_t;

void fmd_vector_init(fmd_vector_t **vec, int_t initial_capacity, fmd_fstruct_t *f);
void fmd_vector_free(fmd_vector_t *vec);
void fmd_vector_grow(fmd_vector_t *vec);
void fmd_vector_shrink(fmd_vector_t *vec);
void fmd_vector_append(fmd_vector_t *vec, void *value);
void fmd_vector_pop(fmd_vector_t *vec, void **item);
void fmd_vector_insert(fmd_vector_t *vec, int_t index, void *value);
void fmd_vector_delete(fmd_vector_t *vec, int_t index, void **item);
void fmd_qs_helper_(void **arr, int_t left, int_t right, fcomp comp_f);
void fmd_vector_sort(fmd_vector_t *vec);
// other functions to wrap in a fstruct
int fmd_vector_comp(fmd_vector_t *v1, fmd_vector_t *v2);
uint_t fmd_vector_hash(fmd_vector_t *v);
fmd_vector_t *fmd_vector_copy(fmd_vector_t *v);

static fmd_fstruct_t fmd_fstruct_vector = {
    fmd_vector_comp,
    fmd_vector_hash,
    fmd_vector_free,
    fmd_vector_copy
};

#endif //FMD_FMD_VECTOR_H
