/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_table.h is part of fmd
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
#ifndef FMD_FMD_TABLE_H
#define FMD_FMD_TABLE_H

#include "fmd_common.h"
#include "fmd_tree.h"

#define FMD_HT_INIT_SIZE 53
#define FMD_REHASH_FACTOR 0.50

typedef struct fmd_table_ {
    int_t size;
    int_t capacity;
    fmd_fstruct_t *key_f;
    fmd_fstruct_t *val_f;
    fmd_tree_node_t **roots;
    int_t *items_per_bucket;
} fmd_table_t;

// internal functions
void fmd_table_rehash_helper_(void* key, void* value, void *ht);
// exposed API
void fmd_table_init(fmd_table_t **table, int_t capacity, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
void fmd_table_free(fmd_table_t *table);
bool fmd_table_insert(fmd_table_t *table, void *key, void *value);
bool fmd_table_lookup(fmd_table_t *table, void *key, void **value);
void fmd_table_rehash(fmd_table_t **table, int_t new_capacity);
int fmd_table_comp(fmd_table_t *t1, fmd_table_t *t2);
fmd_table_t *fmd_table_copy(fmd_table_t *table);
uint_t fmd_table_hash(fmd_table_t *table);
void fmd_table_traverse(fmd_table_t *table, void *p, ftrav_kv f);

static fmd_fstruct_t fmd_fstruct_table = {
        (fcomp) fmd_table_comp,
        (fhash) fmd_table_hash,
        (ffree) fmd_table_free,
        (fcopy) fmd_table_copy
};

#endif //FMD_FMD_TABLE_H
