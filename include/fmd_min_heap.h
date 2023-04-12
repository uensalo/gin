#ifndef FMD_FMD_MIN_HEAP_H
#define FMD_FMD_MIN_HEAP_H

#include "fmd_common.h"

typedef struct fmd_min_heap_kv_s {
    void *key;
    void *value;
} fmd_min_heap_kv_t;

typedef struct fmd_min_heap_s {
    fmd_min_heap_kv_t **data;
    int_t capacity;
    int_t size;
    fmd_fstruct_t *key_f;
    fmd_fstruct_t *val_f;
} fmd_min_heap_t;

void fmd_min_heap_init(fmd_min_heap_t **heap, int_t capacity, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
void fmd_min_heap_free(fmd_min_heap_t *heap);
bool fmd_min_heap_push(fmd_min_heap_t *heap, void *key, void* value);
bool fmd_min_heap_pop(fmd_min_heap_t *heap, void **key, void **value);
bool fmd_min_heap_peek(fmd_min_heap_t *heap, void **key, void **value);

// functions to maintain the heap property
void fmd_min_heap_sift_up(fmd_min_heap_t *heap, int_t index);
void fmd_min_heap_sift_down(fmd_min_heap_t *heap, int_t index);

// other functions to wrap this into an fstruct
fmd_min_heap_t *fmd_min_heap_copy(fmd_min_heap_t *heap);
uint_t fmd_min_heap_hash(fmd_min_heap_t *heap);
int fmd_min_heap_comp(fmd_min_heap_t *h1, fmd_min_heap_t *h2);

fmd_fstruct_t fmd_fstruct_min_heap = {
    fmd_min_heap_comp,
    fmd_min_heap_hash,
    fmd_min_heap_free,
    fmd_min_heap_copy
};

#endif //FMD_FMD_MIN_HEAP_H
