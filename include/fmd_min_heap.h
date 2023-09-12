/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_min_heap.h is part of fmd
 *
 * fmd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * fmd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FMD_FMD_MIN_HEAP_H
#define FMD_FMD_MIN_HEAP_H

#include "fmd_common.h"

typedef struct fmd_min_heap_kv_s {
    void *key;
    void *value;
} fmd_min_heap_kv_t;

typedef struct fmd_min_heap_s {
    fmd_min_heap_kv_t *data;
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

static fmd_fstruct_t fmd_fstruct_min_heap = {
        (fcomp) fmd_min_heap_comp,
        (fhash) fmd_min_heap_hash,
        (ffree) fmd_min_heap_free,
        (fcopy) fmd_min_heap_copy
};

#endif //FMD_FMD_MIN_HEAP_H
