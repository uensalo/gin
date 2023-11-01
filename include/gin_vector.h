/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_vector.h is part of gin
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
#ifndef GIN_GIN_VECTOR_H
#define GIN_GIN_VECTOR_H
#include "gin_common.h"

#define GIN_VECTOR_INIT_SIZE 8

typedef struct gin_vector_ {
    void **data;
    int_t size;
    int_t capacity;
    gin_fstruct_t *f;
} gin_vector_t;

void gin_vector_init(gin_vector_t **vec, int_t initial_capacity, gin_fstruct_t *f);
void gin_vector_free(gin_vector_t *vec);
void gin_vector_free_disown(gin_vector_t *vec);
void gin_vector_grow(gin_vector_t *vec);
void gin_vector_shrink(gin_vector_t *vec);
void gin_vector_fit(gin_vector_t *vec);
void gin_vector_append(gin_vector_t *vec, void *value);
void gin_vector_pop(gin_vector_t *vec, void **item);
void gin_vector_insert(gin_vector_t *vec, int_t index, void *value);
void gin_vector_delete(gin_vector_t *vec, int_t index, void **item);
void gin_qs_helper_(void **arr, int_t left, int_t right, fcomp comp_f);
void gin_qs_arg_helper_(void **arr, void **args, int_t left, int_t right, fcomp comp_f);
void gin_vector_sort(gin_vector_t *vec);
void gin_vector_argsort(gin_vector_t **args, gin_vector_t *vec);
// other functions to wrap in a fstruct
int gin_vector_comp(gin_vector_t *v1, gin_vector_t *v2);
uint_t gin_vector_hash(gin_vector_t *v);
gin_vector_t *gin_vector_copy(gin_vector_t *v);

static gin_fstruct_t gin_fstruct_vector = {
        (fcomp) gin_vector_comp,
        (fhash) gin_vector_hash,
        (ffree) gin_vector_free,
        (fcopy) gin_vector_copy
};

#endif //GIN_GIN_VECTOR_H
