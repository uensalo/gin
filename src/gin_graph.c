/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_graph.c is part of gin
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
#include "gin_graph.h"

void gin_vertex_init(gin_vertex_t **vertex, vid_t id, gin_string_t *label) {
    gin_vertex_t *v = calloc(1, sizeof(gin_vertex_t));
    if(!v) {
        *vertex = NULL;
        return;
    }
    v->id = id;
    v->label = label;
    *vertex = v;
}

void gin_vertex_free(gin_vertex_t *vertex) {
    gin_string_free(vertex->label);
    free(vertex);
}

uint_t gin_vertex_hash(gin_vertex_t *vertex) {
    // djb2
    uint_t hash = 5381;
    hash = ((hash << 5) + hash) + (size_t)vertex->id;
    uint_t label_hash = gin_string_hash(vertex->label);
    hash = ((hash << 5) + hash) + label_hash;
    return hash;
}

int gin_vertex_comp(gin_vertex_t *v1, gin_vertex_t *v2) {
    if (v1->id != v2->id) {
        return v1->id < v2->id ? -1 : 1;
    }
    return gin_string_comp(v1->label, v2->label);
}

gin_vertex_t *gin_vertex_copy(gin_vertex_t *vertex) {
    gin_vertex_t *copy = calloc(1,sizeof(gin_vertex_t));
    if(!copy) return NULL;
    copy->id = vertex->id;
    copy->label = gin_string_copy(vertex->label);
    return copy;
}

void gin_graph_init(gin_graph_t **graph) {
    *graph = calloc(1, sizeof(gin_graph_t));
    gin_vector_init(&(*graph)->vertex_list, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    gin_table_init(&(*graph)->vertices, GIN_HT_INIT_SIZE, &prm_fstruct, &gin_fstruct_vertex);
    gin_table_init(&(*graph)->incoming_neighbors, GIN_HT_INIT_SIZE, &prm_fstruct, &gin_fstruct_vector);
    gin_table_init(&(*graph)->outgoing_neighbors, GIN_HT_INIT_SIZE, &prm_fstruct, &gin_fstruct_vector);
}

void gin_graph_free(gin_graph_t *graph) {
    gin_vector_free(graph->vertex_list);
    gin_table_free(graph->vertices);
    gin_table_free(graph->incoming_neighbors);
    gin_table_free(graph->outgoing_neighbors);
    free(graph);
}

void gin_graph_insert_vertex(gin_graph_t *graph, vid_t id, gin_string_t *label) {
    gin_vertex_t *vertex;
    gin_vertex_init(&vertex, id, label);
    gin_table_insert(graph->vertices, (void*)id, vertex);
    gin_vector_append(graph->vertex_list, vertex);

    gin_vector_t *incoming_neighbors;
    gin_vector_init(&incoming_neighbors, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    gin_table_insert(graph->incoming_neighbors, (void*)id, incoming_neighbors);

    gin_vector_t *outgoing_neighbors;
    gin_vector_init(&outgoing_neighbors, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    gin_table_insert(graph->outgoing_neighbors, (void*)id, outgoing_neighbors);
}

void gin_graph_insert_edge(gin_graph_t *graph, vid_t source, vid_t destination) {
    void *incoming_neighbors;
    gin_table_lookup(graph->incoming_neighbors, (void*)destination, &incoming_neighbors);
    gin_vector_append((gin_vector_t*)incoming_neighbors, (void*)source);

    void *outgoing_neigbors;
    gin_table_lookup(graph->outgoing_neighbors, (void*)source, &outgoing_neigbors);
    gin_vector_append((gin_vector_t*)outgoing_neigbors, (void*)destination);

    ++graph->no_edges;
}

uint_t gin_graph_hash(gin_graph_t *graph) {
    uint_t vertices_hash = gin_table_hash(graph->vertices);
    uint_t edges_hash = gin_table_hash(graph->incoming_neighbors);
    return (vertices_hash * 31) + edges_hash;
}

int gin_graph_comp(gin_graph_t *g1, gin_graph_t *g2) {
    int vertices_comp = gin_table_comp(g1->vertices, g2->vertices);
    if (vertices_comp != 0) {
        return vertices_comp;
    }
    return gin_table_comp(g1->incoming_neighbors, g2->incoming_neighbors);
}

gin_graph_t *gin_graph_copy(gin_graph_t *graph) {
    gin_graph_t *copy = calloc(1, sizeof(gin_graph_t));
    copy->vertex_list = gin_vector_copy(graph->vertex_list);
    copy->vertices = gin_table_copy(graph->vertices);
    copy->incoming_neighbors = gin_table_copy(graph->incoming_neighbors);
    copy->outgoing_neighbors = gin_table_copy(graph->outgoing_neighbors);
    return copy;
}

int gin_kmer_kv_comp(gin_kmer_kv_t *k1, gin_kmer_kv_t *k2) {
    int strcmpval = gin_string_comp(k1->str, k2->str);
    return strcmpval ? strcmpval : (int)(k1->vid - k2->vid ?  k1->vid - k2->vid : k1->offset - k2->offset);
}

uint_t gin_kmer_kv_hash(gin_kmer_kv_t *k) {
    return prm_hash_f((void*)k->vid) ^ prm_hash_f((void*)k->offset) ^ gin_string_hash(k->str);
}

void gin_kmer_kv_free(gin_kmer_kv_t *k) {
    if(k) {
        free(k);
    }
}

gin_kmer_kv_t *gin_kmer_kv_copy(gin_kmer_kv_t *k) {
    gin_kmer_kv_t *retval = NULL;
    if(k) {
        retval = calloc(1, sizeof(gin_kmer_kv_t));
        retval->str = gin_string_copy(k->str);
        retval->vid = k->vid;
        retval->offset = k->offset;
    }
    return retval;
}

void gin_graph_kmer_locations_helper_trav(void *key, void *value, void *params) {
    gin_vector_t *v = (gin_vector_t*)params;
    gin_vector_t *k = value;
    for(int_t i = 0; i < k->size; i++) {
        gin_kmer_kv_t *kv = k->data[i];
        kv->str = key;
        gin_vector_append(v, kv);
    }
}

void gin_graph_kmer_locations(gin_graph_t *graph, int_t k, gin_vector_t **kmers, gin_table_t **kmer_table) {
    int_t V = graph->vertex_list->size;
    gin_table_t *set;
    gin_table_init(&set, GIN_HT_INIT_SIZE, &gin_fstruct_string, &gin_fstruct_vector);
    for(int_t i = 0; i < V; i++) {
        gin_vertex_t *vertex = graph->vertex_list->data[i];
        for(int_t j = 0; j < vertex->label->size - k + 1; j++) {
            gin_string_t *kmer;
            gin_string_substring(vertex->label, j, j + k, &kmer);
            gin_kmer_kv_t *loc = calloc(1, sizeof(gin_kmer_kv_t));
            loc->vid = i;
            loc->offset = j;
            gin_vector_t *locs;
            void *_;
            bool exists = gin_table_lookup(set, kmer, &_);
            if(!exists) {
                gin_vector_init(&locs, GIN_VECTOR_INIT_SIZE, &gin_fstruct_kmer_kv);
                gin_vector_append(locs, loc);
                gin_string_t *set_kmer = gin_string_copy(kmer);
                gin_table_insert(set, set_kmer, locs);
            } else {
                locs = (gin_vector_t*)_;
                gin_vector_append(locs, loc);
            }
            gin_string_free(kmer);
        }

        // check the outgoing neighbors
        gin_vector_t *outgoing;
        gin_table_lookup(graph->outgoing_neighbors, (void*)vertex->id, &outgoing);
        typedef struct kmer_dfs_ {
            gin_vertex_t *v;
            gin_string_t *str;
            int_t remaining;
            vid_t source_vid;
            int_t source_offset;
        } kmer_dfs_t;
        gin_vector_t *stack;
        gin_vector_init(&stack, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
        // seed the stack
        for(int_t j = 0; j < outgoing->size; j++) {
            gin_vertex_t *nvertex = graph->vertex_list->data[(int_t)outgoing->data[j]];
            int_t bound = MIN2(k-1, vertex->label->size);
            for(int_t q = 1; q <= bound; q++) {
                gin_string_t *suf;
                gin_string_substring(vertex->label, vertex->label->size-q, vertex->label->size, &suf);
                int_t remaining = k - q;
                kmer_dfs_t *rec = calloc(1, sizeof(kmer_dfs_t));
                rec->remaining = remaining;
                rec->str = suf;
                rec->v = nvertex;
                rec->source_vid = i;
                rec->source_offset = vertex->label->size-q;
                gin_vector_append(stack, rec);
            }
        }

        while(stack->size) {
            void *_;
            gin_vector_pop(stack, &_);
            kmer_dfs_t *rec = _;
            int_t n_traversable = MIN2(rec->v->label->size, rec->remaining);
            gin_string_t *prefix;
            gin_string_substring(rec->v->label, 0, n_traversable, &prefix);
            gin_string_concat_mut(rec->str, prefix);
            gin_string_free(prefix);
            if(n_traversable == rec->remaining) {
                gin_kmer_kv_t *loc = calloc(1,sizeof(gin_kmer_kv_t));
                loc->vid = rec->source_vid;
                loc->offset = rec->source_offset;
                bool exists = gin_table_lookup(set, rec->str, &_);
                gin_vector_t *locs;
                if(exists) {
                    locs = _;
                    gin_vector_append(locs,loc);
                } else {
                    gin_vector_init(&locs, GIN_VECTOR_INIT_SIZE, &gin_fstruct_kmer_kv);
                    gin_vector_append(locs, loc);
                    gin_string_t *set_kmer = gin_string_copy(rec->str);
                    if(!set_kmer->seq[0]) {
                        printf("wtf\n");
                    }
                    gin_table_insert(set, set_kmer, locs);
                }
                gin_string_free(rec->str);
                free(rec);
            } else {
                int_t remaining = rec->remaining - n_traversable;
                gin_vector_t *outgoingv;
                gin_table_lookup(graph->outgoing_neighbors, (void*)rec->v->id, &_);
                outgoingv = _;
                if(outgoingv->size) {
                    for (int_t j = 0; j < outgoingv->size; j++) {
                        gin_vertex_t *nvertex = graph->vertex_list->data[(int_t)outgoingv->data[j]];
                        kmer_dfs_t *nrec = calloc(1, sizeof(kmer_dfs_t));
                        nrec->v = nvertex;
                        nrec->str = gin_string_copy(rec->str);
                        nrec->remaining = remaining;
                        nrec->source_vid = rec->source_vid;
                        nrec->source_offset = rec->source_offset;
                        gin_vector_append(stack, nrec);
                    }
                } else { // no outgoing vertices
                    gin_string_free(rec->str);
                    free(rec);
                }
            }
        }
        gin_vector_free(stack);
    }
    // gather all strings into a vector and sort
    gin_vector_t *kmer_vec;
    gin_vector_init(&kmer_vec, set->size, &gin_fstruct_kmer_kv);
    gin_table_traverse(set, kmer_vec, gin_graph_kmer_locations_helper_trav);
    *kmers = kmer_vec;
    *kmer_table = set;
}

void gin_graph_find(gin_graph_t *graph, gin_string_t *string, gin_vector_t **origins) {
    int_t V = graph->vertex_list->size;
    // "compile" the string pattern once
    int_t *lps;
    gin_string_kmp_lps(string, &lps);
    gin_vector_t **vertex_hits = calloc(V, sizeof(gin_vector_t*));
    // first, compute matches without vertex extensions
    #pragma omp parallel for default(none) shared(V,graph,vertex_hits,prm_fstruct,lps,string)
    for(int_t i = 0; i < V; i++) {
        gin_vertex_t *vertex = graph->vertex_list->data[i];
        gin_vector_init(&vertex_hits[i], GIN_VECTOR_INIT_SIZE, &prm_fstruct);
        if(vertex->label->size >= string->size) {
            int_t *offsets;
            int_t no_offsets;
            gin_string_kmp_search(vertex->label, string, lps, &offsets, &no_offsets);
            for (int_t j = 0; j < no_offsets; j++) {
                gin_vector_append(vertex_hits[i], (void*)offsets[j]);
            }
            free(offsets);
        }
    }
    free(lps);
    // then, compute matches only with vertex extensions
    typedef struct graph_find_dfs_ {
        gin_vertex_t *v;
        int_t pos;
    } graph_find_dfs_t;
    #pragma omp parallel for default(none) shared(V,graph,vertex_hits,prm_fstruct,string)
    for(int_t i = 0; i < V; i++) {
        gin_vertex_t *vertex = graph->vertex_list->data[i];
        // find the location of a match of the prefix of the pattern to a suffix of the label
        // the suffix of the label is shorter by construction
        int_t query_prefix_match_length = 0;
        for (int_t j = MAX2(0, vertex->label->size - string->size + 1); j < vertex->label->size; j++) {
            int_t m = 0;
            while (j + m < vertex->label->size && m < string->size && vertex->label->seq[j + m] == string->seq[m]) {
                m++;
            }
            if (j + m == vertex->label->size) {
                query_prefix_match_length = m;
                break;
            }
        }
        if(!query_prefix_match_length) continue;
        int_t pos = query_prefix_match_length; // position on the current match
        gin_vector_t *stack;
        gin_vector_init(&stack, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
        // add all the outgoing neighbors;
        gin_vector_t *n;
        gin_table_lookup(graph->outgoing_neighbors, (void*)vertex->id, &n);
        for(int_t j = 0; j < n->size; j++) {
            graph_find_dfs_t *rec = malloc(sizeof(graph_find_dfs_t));
            rec->v = graph->vertex_list->data[(int_t)n->data[j]];
            rec->pos = pos;
            gin_vector_append(stack, rec);
        }
        while(stack->size) {
            void* item;
            gin_vector_pop(stack, &item);
            graph_find_dfs_t *rec = item;
            int_t n_traversable = MIN2(rec->v->label->size, string->size-rec->pos);
            bool matched = true;
            for(int_t j = 0; j < n_traversable; j++) {
                if(string->seq[j+rec->pos] != rec->v->label->seq[j]) {
                    free(rec);
                    matched = false;
                    break;
                }
            }
            if(matched) {
                if(rec->pos + n_traversable == string->size) { // match exhausted
                    gin_vector_append(vertex_hits[i], (void*)(vertex->label->size - query_prefix_match_length));
                    free(rec);
                } else { // match needs to be extended
                    gin_table_lookup(graph->outgoing_neighbors, (void*)rec->v->id, &n);
                    for(int_t j = 0; j < n->size; j++) {
                        graph_find_dfs_t *nrec = malloc(sizeof(graph_find_dfs_t));
                        nrec->v = graph->vertex_list->data[(int_t)n->data[j]];
                        nrec->pos = rec->pos + n_traversable;
                        gin_vector_append(stack, nrec);
                    }
                }
            }
        }
        gin_vector_free(stack);
    }

    gin_vector_init(origins, GIN_VECTOR_INIT_SIZE, &gin_fstruct_kmer_kv);
    // cleanup and return

    #pragma omp parallel for default(none) shared(V,vertex_hits)
    for(int_t i = 0; i < V; i++) {
        gin_vector_sort(vertex_hits[i]);
    }

    for(int_t i = 0; i < V; i++) {
        for(int_t j = 0; j < vertex_hits[i]->size; j++) {
            gin_kmer_kv_t *kv = malloc(sizeof(gin_kmer_kv_t));
            kv->str = string;
            kv->vid = i;
            kv->offset = (int_t)vertex_hits[i]->data[j];
            gin_vector_append(*origins, kv);
        }
        gin_vector_free(vertex_hits[i]);
    }
    free(vertex_hits);
}