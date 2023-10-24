/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_constraint_set.h is part of gin
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
#ifndef GIN_GIN_CONSTRAINT_SET_H
#define GIN_GIN_CONSTRAINT_SET_H

#include "gin_common.h"
#include "gin_graph.h"

typedef struct gin_constraint_set_ {
    gin_string_t *str;
    gin_vector_t *vertices;
} gin_constraint_set_t;

void gin_constraint_set_enumerate(gin_vector_t **constraint_sets, gin_graph_t *graph, int_t depth, bool multiple_vertex_span);
void gin_constraint_set_free(gin_constraint_set_t *cs);
uint_t gin_constraint_set_hash(gin_constraint_set_t *c1);
int gin_constraint_set_comp(gin_constraint_set_t *c1, gin_constraint_set_t *c2);
void *gin_constraint_set_copy(gin_constraint_set_t *c1);

static gin_fstruct_t gin_fstruct_cs = {
        (fcomp) gin_constraint_set_comp,
        (fhash) gin_constraint_set_hash,
        (ffree) gin_constraint_set_free,
        (fcopy) gin_constraint_set_copy
};

typedef struct gin_graph_flatten_alphabet_ {
    gin_vector_t *alphabet;
    gin_table_t *idx2char;
    gin_table_t *char2idx;
} gin_graph_flatten_alphabet_t;
void gin_constraint_set_flatten_alphabet(void *key, void *value, void *params);
typedef struct gin_graph_ecs_ {
    int_t pos;
    vid_t head_vid;
    vid_t end_vid;
} gin_graph_ecs_t;
gin_table_t *gin_constraint_set_extract(gin_graph_t *graph, int_t max_depth, bool multiple_vertex_span);
void gin_constraint_set_extract_helper(gin_vector_t *paths,
                                              gin_string_t *prefix,
                                              gin_table_t *constraint_sets,
                                              gin_vector_t *alphabet,
                                              gin_table_t *char2idx,
                                              gin_table_t *idx2char,
                                              gin_graph_t *graph,
                                              int_t max_depth,
                                              bool multiple_vertex_span);

void gin_constraint_set_flatten(void *key, void *value, void *p);

void gin_constraint_set_to_vector_helper(void *key, void* value, void *p);
gin_vector_t *gin_constraint_set_to_vector(gin_table_t *constraint_sets);

#endif //GIN_GIN_CONSTRAINT_SET_H
