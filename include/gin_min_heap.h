/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_min_heap.h is part of gin
 *
 * gin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GIN_GIN_MIN_HEAP_H
#define GIN_GIN_MIN_HEAP_H

#include "gin_common.h"

typedef struct gin_min_heap_kv_s {
    void *key;
    void *value;
} gin_min_heap_kv_t;

typedef struct gin_min_heap_s {
    gin_min_heap_kv_t *data;
    int_t capacity;
    int_t size;
    gin_fstruct_t *key_f;
    gin_fstruct_t *val_f;
} gin_min_heap_t;

void gin_min_heap_init(gin_min_heap_t **heap, int_t capacity, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
void gin_min_heap_free(gin_min_heap_t *heap);
bool gin_min_heap_push(gin_min_heap_t *heap, void *key, void* value);
bool gin_min_heap_pop(gin_min_heap_t *heap, void **key, void **value);
bool gin_min_heap_peek(gin_min_heap_t *heap, void **key, void **value);

// functions to maintain the heap property
void gin_min_heap_sift_up(gin_min_heap_t *heap, int_t index);
void gin_min_heap_sift_down(gin_min_heap_t *heap, int_t index);

// other functions to wrap this into an fstruct
gin_min_heap_t *gin_min_heap_copy(gin_min_heap_t *heap);
uint_t gin_min_heap_hash(gin_min_heap_t *heap);
int gin_min_heap_comp(gin_min_heap_t *h1, gin_min_heap_t *h2);

static gin_fstruct_t gin_fstruct_min_heap = {
        (fcomp) gin_min_heap_comp,
        (fhash) gin_min_heap_hash,
        (ffree) gin_min_heap_free,
        (fcopy) gin_min_heap_copy
};

#endif //GIN_GIN_MIN_HEAP_H
