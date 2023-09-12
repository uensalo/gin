/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_constraint_set.h is part of fmd
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
#ifndef FMD_FMD_CONSTRAINT_SET_H
#define FMD_FMD_CONSTRAINT_SET_H

#include "fmd_common.h"
#include "fmd_graph.h"

typedef struct fmd_constraint_set_ {
    fmd_string_t *str;
    fmd_vector_t *vertices;
} fmd_constraint_set_t;

void fmd_constraint_set_enumerate(fmd_vector_t **constraint_sets, fmd_graph_t *graph, int_t depth, bool multiple_vertex_span);
void fmd_constraint_set_free(fmd_constraint_set_t *cs);
uint_t fmd_constraint_set_hash(fmd_constraint_set_t *c1);
int fmd_constraint_set_comp(fmd_constraint_set_t *c1, fmd_constraint_set_t *c2);
void *fmd_constraint_set_copy(fmd_constraint_set_t *c1);

static fmd_fstruct_t fmd_fstruct_cs = {
        (fcomp) fmd_constraint_set_comp,
        (fhash) fmd_constraint_set_hash,
        (ffree) fmd_constraint_set_free,
        (fcopy) fmd_constraint_set_copy
};

typedef struct fmd_graph_flatten_alphabet_ {
    fmd_vector_t *alphabet;
    fmd_table_t *idx2char;
    fmd_table_t *char2idx;
} fmd_graph_flatten_alphabet_t;
void fmd_constraint_set_flatten_alphabet(void *key, void *value, void *params);
typedef struct fmd_graph_ecs_ {
    int_t pos;
    vid_t head_vid;
    vid_t end_vid;
} fmd_graph_ecs_t;
fmd_table_t *fmd_constraint_set_extract(fmd_graph_t *graph, int_t max_depth, bool multiple_vertex_span);
void fmd_constraint_set_extract_helper(fmd_vector_t *paths,
                                              fmd_string_t *prefix,
                                              fmd_table_t *constraint_sets,
                                              fmd_vector_t *alphabet,
                                              fmd_table_t *char2idx,
                                              fmd_table_t *idx2char,
                                              fmd_graph_t *graph,
                                              int_t max_depth,
                                              bool multiple_vertex_span);

void fmd_constraint_set_flatten(void *key, void *value, void *p);

void fmd_constraint_set_to_vector_helper(void *key, void* value, void *p);
fmd_vector_t *fmd_constraint_set_to_vector(fmd_table_t *constraint_sets);

#endif //FMD_FMD_CONSTRAINT_SET_H
