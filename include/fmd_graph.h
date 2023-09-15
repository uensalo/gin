/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_graph.h is part of fmd
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
#ifndef FMD_FMD_GRAPH_H
#define FMD_FMD_GRAPH_H

#include <fmd_string.h>
#include <fmd_table.h>
#include <fmd_tree.h>
#include <fmd_vector.h>

typedef struct {
    vid_t id;
    fmd_string_t *label;
} fmd_vertex_t;
void fmd_vertex_init(fmd_vertex_t **vertex, vid_t id, fmd_string_t *label);
void fmd_vertex_free(fmd_vertex_t *vertex);
uint_t fmd_vertex_hash(fmd_vertex_t *vertex);
int fmd_vertex_comp(fmd_vertex_t *v1, fmd_vertex_t *v2);
fmd_vertex_t *fmd_vertex_copy(fmd_vertex_t *vertex);

static fmd_fstruct_t fmd_fstruct_vertex = {
        (fcomp) fmd_vertex_comp,
        (fhash) fmd_vertex_hash,
        (ffree) fmd_vertex_free,
        (fcopy) fmd_vertex_copy
};

typedef struct {
    fmd_vector_t *vertex_list;
    fmd_table_t *vertices;
    fmd_table_t *incoming_neighbors;
    fmd_table_t *outgoing_neighbors;
    int_t no_edges;
} fmd_graph_t;

// Graph functions
void fmd_graph_init(fmd_graph_t **graph);
void fmd_graph_free(fmd_graph_t *graph);
void fmd_graph_insert_vertex(fmd_graph_t *graph, vid_t id, fmd_string_t *label);
void fmd_graph_insert_edge(fmd_graph_t *graph, vid_t source, vid_t destination);
uint_t fmd_graph_hash(fmd_graph_t *graph);
int fmd_graph_comp(fmd_graph_t *g1, fmd_graph_t *g2);
fmd_graph_t *fmd_graph_copy(fmd_graph_t *graph);
static fmd_fstruct_t fmd_fstruct_graph = {
        (fcomp) fmd_graph_comp,
        (fhash) fmd_graph_hash,
        (ffree) fmd_graph_free,
        (fcopy) fmd_graph_copy
};

typedef struct fmd_kmer_kv_ {
    fmd_string_t *str;
    vid_t vid;
    int_t offset;
} fmd_kmer_kv_t;
int fmd_kmer_kv_comp(fmd_kmer_kv_t *k1, fmd_kmer_kv_t *k2);
uint_t fmd_kmer_kv_hash(fmd_kmer_kv_t *k);
void fmd_kmer_kv_free(fmd_kmer_kv_t *k);
fmd_kmer_kv_t *fmd_kmer_kv_copy(fmd_kmer_kv_t *k);

void fmd_graph_kmer_locations_helper_trav(void *key, void *value, void *params);
void fmd_graph_kmer_locations(fmd_graph_t *graph, int_t k, fmd_vector_t **kmers, fmd_table_t **table);

static fmd_fstruct_t fmd_fstruct_kmer_kv = {
        (fcomp) fmd_kmer_kv_comp,
        (fhash) fmd_kmer_kv_hash,
        (ffree) fmd_kmer_kv_free,
        (fcopy) fmd_kmer_kv_copy
};

void fmd_graph_find(fmd_graph_t *graph, fmd_string_t *string, fmd_vector_t **origins);

#endif //FMD_FMD_GRAPH_H
