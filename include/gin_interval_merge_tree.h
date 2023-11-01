/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_interval_merge_tree.h is part of gin
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
#ifndef GIN_GIN_INTERVAL_MERGE_TREE_H
#define GIN_GIN_INTERVAL_MERGE_TREE_H
#include "gin_vector.h"
#include "gin_min_heap.h"

typedef struct gin_imt_interval_{
    int_t lo;
    int_t hi;
} gin_imt_interval_t;

void gin_imt_interval_init(gin_imt_interval_t **i, int_t lo, int_t hi);
void gin_imt_interval_free(gin_imt_interval_t *i);
gin_imt_interval_t *gin_imt_interval_copy(gin_imt_interval_t *i);
uint_t gin_imt_interval_hash(gin_imt_interval_t *i);
inline int gin_imt_interval_comp(gin_imt_interval_t *i1, gin_imt_interval_t *i2); // inline does make a difference :)

static gin_fstruct_t gin_fstruct_imt_interval = {
        (fcomp) gin_imt_interval_comp,
        (fhash) gin_imt_interval_hash,
        (ffree) gin_imt_interval_free,
        (fcopy) gin_imt_interval_copy
};

typedef struct gin_imt_node_ {
    int_t lo;
    int_t hi;
    gin_vector_t *intervals;
    struct gin_imt_node_t *left;
    struct gin_imt_node_t *right;
} gin_imt_node_t;

void gin_imt_node_init(gin_imt_node_t **i, int_t lo, int_t hi, gin_vector_t *intervals, gin_imt_node_t *left, gin_imt_node_t *right);
void gin_imt_node_free(gin_imt_node_t *i);
gin_imt_node_t *gin_imt_node_copy(gin_imt_node_t *i);
uint_t gin_imt_node_hash(gin_imt_node_t *i);
int gin_imt_node_comp(gin_imt_node_t *i1, gin_imt_node_t *i2);

static gin_fstruct_t gin_fstruct_imt_node = {
        (fcomp) gin_imt_node_comp,
        (fhash) gin_imt_node_hash,
        (ffree) gin_imt_node_free,
        (fcopy) gin_imt_node_copy
};

typedef struct gin_imt_ {
    int_t no_keys;
    gin_imt_node_t *root;
} gin_imt_t;

void gin_imt_init(gin_imt_t **i, int_t no_keys, gin_vector_t *kv_interval_pairs);
void gin_imt_free(gin_imt_t *i);
gin_imt_t *gin_imt_copy(gin_imt_t *i);
uint_t gin_imt_hash(gin_imt_t *i);
int gin_imt_comp(gin_imt_t *i1, gin_imt_t *i2);

void gin_imt_query(gin_imt_t *i, int_t start, int_t end, int_t no_max_intervals, gin_vector_t **intervals);
void gin_imt_query_helper(gin_imt_node_t *node, int_t lo, int_t hi, int_t no_max_intervals, int_t *no_cur_intervals, gin_vector_t *merge_list);

void gin_imt_query_legacy(gin_imt_t *i, int_t start, int_t end, gin_vector_t **intervals);
gin_vector_t *gin_imt_query_helper_legacy(gin_imt_node_t *node, int_t lo, int_t hi);

void gin_imt_init_compact_intervals(gin_vector_t *intervals, gin_vector_t **compacted);
gin_imt_node_t *gin_imt_init_helper(int_t lo, int_t hi, gin_vector_t *kv_interval_pairs);
void gin_imt_free_helper(gin_imt_node_t *node);
gin_vector_t *gin_imt_merge_intervals(gin_vector_t *i1, gin_vector_t *i2);
gin_vector_t *gin_imt_multiway_merge_intervals(gin_vector_t *list_of_intervals, int_t no_max_intervals);

gin_imt_node_t *gin_imt_copy_helper(gin_imt_node_t *i);
uint_t gin_imt_hash_helper(gin_imt_node_t *i);
int gin_imt_comp_helper(gin_imt_node_t *i1, gin_imt_node_t *i2);

static gin_fstruct_t gin_fstruct_imt = {
        (fcomp) gin_imt_comp,
        (fhash) gin_imt_hash,
        (ffree) gin_imt_free,
        (fcopy) gin_imt_copy
};

#endif //GIN_GIN_INTERVAL_MERGE_TREE_H
