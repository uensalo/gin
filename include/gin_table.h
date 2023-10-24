/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_table.h is part of gin
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
#ifndef GIN_GIN_TABLE_H
#define GIN_GIN_TABLE_H

#include "gin_common.h"
#include "gin_tree.h"

#define GIN_HT_INIT_SIZE 53
#define GIN_REHASH_FACTOR 0.50

typedef struct gin_table_ {
    int_t size;
    int_t capacity;
    gin_fstruct_t *key_f;
    gin_fstruct_t *val_f;
    gin_tree_node_t **roots;
    int_t *items_per_bucket;
} gin_table_t;

// internal functions
void gin_table_rehash_helper_(void* key, void* value, void *ht);
// exposed API
void gin_table_init(gin_table_t **table, int_t capacity, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
void gin_table_free(gin_table_t *table);
bool gin_table_insert(gin_table_t *table, void *key, void *value);
bool gin_table_lookup(gin_table_t *table, void *key, void **value);
void gin_table_rehash(gin_table_t **table, int_t new_capacity);
int gin_table_comp(gin_table_t *t1, gin_table_t *t2);
gin_table_t *gin_table_copy(gin_table_t *table);
uint_t gin_table_hash(gin_table_t *table);
void gin_table_traverse(gin_table_t *table, void *p, ftrav_kv f);

static gin_fstruct_t gin_fstruct_table = {
        (fcomp) gin_table_comp,
        (fhash) gin_table_hash,
        (ffree) gin_table_free,
        (fcopy) gin_table_copy
};

#endif //GIN_GIN_TABLE_H
