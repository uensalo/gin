/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_vector.h is part of fmd
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
void fmd_vector_free_disown(fmd_vector_t *vec);
void fmd_vector_grow(fmd_vector_t *vec);
void fmd_vector_shrink(fmd_vector_t *vec);
void fmd_vector_fit(fmd_vector_t *vec);
void fmd_vector_append(fmd_vector_t *vec, void *value);
void fmd_vector_pop(fmd_vector_t *vec, void **item);
void fmd_vector_insert(fmd_vector_t *vec, int_t index, void *value);
void fmd_vector_delete(fmd_vector_t *vec, int_t index, void **item);
void fmd_qs_helper_(void **arr, int_t left, int_t right, fcomp comp_f);
void fmd_qs_arg_helper_(void **arr, void **args, int_t left, int_t right, fcomp comp_f);
void fmd_vector_sort(fmd_vector_t *vec);
void fmd_vector_argsort(fmd_vector_t **args, fmd_vector_t *vec);
// other functions to wrap in a fstruct
int fmd_vector_comp(fmd_vector_t *v1, fmd_vector_t *v2);
uint_t fmd_vector_hash(fmd_vector_t *v);
fmd_vector_t *fmd_vector_copy(fmd_vector_t *v);

static fmd_fstruct_t fmd_fstruct_vector = {
        (fcomp) fmd_vector_comp,
        (fhash) fmd_vector_hash,
        (ffree) fmd_vector_free,
        (fcopy) fmd_vector_copy
};

#endif //FMD_FMD_VECTOR_H
