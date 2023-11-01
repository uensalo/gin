/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_graph.h is part of gin
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
#ifndef GIN_GIN_GRAPH_H
#define GIN_GIN_GRAPH_H

#include <gin_string.h>
#include <gin_table.h>
#include <gin_tree.h>
#include <gin_vector.h>

typedef struct {
    vid_t id;
    gin_string_t *label;
} gin_vertex_t;
void gin_vertex_init(gin_vertex_t **vertex, vid_t id, gin_string_t *label);
void gin_vertex_free(gin_vertex_t *vertex);
uint_t gin_vertex_hash(gin_vertex_t *vertex);
int gin_vertex_comp(gin_vertex_t *v1, gin_vertex_t *v2);
gin_vertex_t *gin_vertex_copy(gin_vertex_t *vertex);

static gin_fstruct_t gin_fstruct_vertex = {
        (fcomp) gin_vertex_comp,
        (fhash) gin_vertex_hash,
        (ffree) gin_vertex_free,
        (fcopy) gin_vertex_copy
};

typedef struct {
    gin_vector_t *vertex_list;
    gin_table_t *vertices;
    gin_table_t *incoming_neighbors;
    gin_table_t *outgoing_neighbors;
    int_t no_edges;
} gin_graph_t;

// Graph functions
void gin_graph_init(gin_graph_t **graph);
void gin_graph_free(gin_graph_t *graph);
void gin_graph_insert_vertex(gin_graph_t *graph, vid_t id, gin_string_t *label);
void gin_graph_insert_edge(gin_graph_t *graph, vid_t source, vid_t destination);
uint_t gin_graph_hash(gin_graph_t *graph);
int gin_graph_comp(gin_graph_t *g1, gin_graph_t *g2);
gin_graph_t *gin_graph_copy(gin_graph_t *graph);
static gin_fstruct_t gin_fstruct_graph = {
        (fcomp) gin_graph_comp,
        (fhash) gin_graph_hash,
        (ffree) gin_graph_free,
        (fcopy) gin_graph_copy
};

typedef struct gin_kmer_kv_ {
    gin_string_t *str;
    vid_t vid;
    int_t offset;
} gin_kmer_kv_t;
int gin_kmer_kv_comp(gin_kmer_kv_t *k1, gin_kmer_kv_t *k2);
uint_t gin_kmer_kv_hash(gin_kmer_kv_t *k);
void gin_kmer_kv_free(gin_kmer_kv_t *k);
gin_kmer_kv_t *gin_kmer_kv_copy(gin_kmer_kv_t *k);

void gin_graph_kmer_locations_helper_trav(void *key, void *value, void *params);
void gin_graph_kmer_locations(gin_graph_t *graph, int_t k, gin_vector_t **kmers, gin_table_t **table);

static gin_fstruct_t gin_fstruct_kmer_kv = {
        (fcomp) gin_kmer_kv_comp,
        (fhash) gin_kmer_kv_hash,
        (ffree) gin_kmer_kv_free,
        (fcopy) gin_kmer_kv_copy
};

void gin_graph_find(gin_graph_t *graph, gin_string_t *string, gin_vector_t **origins);

#endif //GIN_GIN_GRAPH_H
