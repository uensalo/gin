#include "fmd_graph.h"

void fmd_vertex_init(fmd_vertex_t **vertex, vid_t id, fmd_string_t *label) {
    fmd_vertex_t *v = calloc(1, sizeof(fmd_vertex_t));
    if(!v) {
        *vertex = NULL;
        return;
    }
    v->id = id;
    v->label = label;
    *vertex = v;
}

void fmd_vertex_free(fmd_vertex_t *vertex) {
    fmd_string_free(vertex->label);
    free(vertex);
}

uint_t fmd_vertex_hash(fmd_vertex_t *vertex) {
    // djb2
    uint_t hash = 5381;
    hash = ((hash << 5) + hash) + (size_t)vertex->id;
    uint_t label_hash = fmd_string_hash(vertex->label);
    hash = ((hash << 5) + hash) + label_hash;
    return hash;
}

int fmd_vertex_comp(fmd_vertex_t *v1, fmd_vertex_t *v2) {
    if (v1->id != v2->id) {
        return v1->id < v2->id ? -1 : 1;
    }
    return fmd_string_comp(v1->label, v2->label);
}

fmd_vertex_t *fmd_vertex_copy(fmd_vertex_t *vertex) {
    fmd_vertex_t *copy = calloc(1,sizeof(fmd_vertex_t));
    if(!copy) return NULL;
    copy->id = vertex->id;
    copy->label = fmd_string_copy(vertex->label);
    return copy;
}

void fmd_graph_init(fmd_graph_t **graph) {
    *graph = calloc(1, sizeof(fmd_graph_t));
    fmd_vector_init(&(*graph)->vertex_list, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_init(&(*graph)->vertices, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vertex);
    fmd_table_init(&(*graph)->incoming_neighbors, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vector);
    fmd_table_init(&(*graph)->outgoing_neighbors, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vector);
}

void fmd_graph_free(fmd_graph_t *graph) {
    fmd_vector_free(graph->vertex_list);
    fmd_table_free(graph->vertices);
    fmd_table_free(graph->incoming_neighbors);
    fmd_table_free(graph->outgoing_neighbors);
    free(graph);
}

void fmd_graph_insert_vertex(fmd_graph_t *graph, vid_t id, fmd_string_t *label) {
    fmd_vertex_t *vertex;
    fmd_vertex_init(&vertex, id, label);
    fmd_table_insert(graph->vertices, (void*)id, vertex);
    fmd_vector_append(graph->vertex_list, vertex);

    fmd_vector_t *incoming_neighbors;
    fmd_vector_init(&incoming_neighbors, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_insert(graph->incoming_neighbors, (void*)id, incoming_neighbors);

    fmd_vector_t *outgoing_neighbors;
    fmd_vector_init(&outgoing_neighbors, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_insert(graph->outgoing_neighbors, (void*)id, outgoing_neighbors);
}

void fmd_graph_insert_edge(fmd_graph_t *graph, vid_t source, vid_t destination) {
    void *incoming_neighbors;
    fmd_table_lookup(graph->incoming_neighbors, (void*)destination, &incoming_neighbors);
    fmd_vector_append((fmd_vector_t*)incoming_neighbors, (void*)source);

    void *outgoing_neigbors;
    fmd_table_lookup(graph->outgoing_neighbors, (void*)source, &outgoing_neigbors);
    fmd_vector_append((fmd_vector_t*)outgoing_neigbors, (void*)destination);

    ++graph->no_edges;
}

uint_t fmd_graph_hash(fmd_graph_t *graph) {
    uint_t vertices_hash = fmd_table_hash(graph->vertices);
    uint_t edges_hash = fmd_table_hash(graph->incoming_neighbors);
    return (vertices_hash * 31) + edges_hash;
}

int fmd_graph_comp(fmd_graph_t *g1, fmd_graph_t *g2) {
    int vertices_comp = fmd_table_comp(g1->vertices, g2->vertices);
    if (vertices_comp != 0) {
        return vertices_comp;
    }
    return fmd_table_comp(g1->incoming_neighbors, g2->incoming_neighbors);
}

fmd_graph_t *fmd_graph_copy(fmd_graph_t *graph) {
    fmd_graph_t *copy = calloc(1, sizeof(fmd_graph_t));
    copy->vertex_list = fmd_vector_copy(graph->vertex_list);
    copy->vertices = fmd_table_copy(graph->vertices);
    copy->incoming_neighbors = fmd_table_copy(graph->incoming_neighbors);
    copy->outgoing_neighbors = fmd_table_copy(graph->outgoing_neighbors);
    return copy;
}

int fmd_kmer_kv_comp(fmd_kmer_kv_t *k1, fmd_kmer_kv_t *k2) {
    int strcmpval = fmd_string_comp(k1->str, k2->str);
    return strcmpval ? strcmpval : (int)(k1->vid - k2->vid ?  k1->vid - k2->vid : k1->offset - k2->offset);
}

uint_t fmd_kmer_kv_hash(fmd_kmer_kv_t *k) {
    return prm_hash_f((void*)k->vid) ^ prm_hash_f((void*)k->offset) ^ fmd_string_hash(k->str);
}

void fmd_kmer_kv_free(fmd_kmer_kv_t *k) {
    if(k) {
        free(k);
    }
}

fmd_kmer_kv_t *fmd_kmer_kv_copy(fmd_kmer_kv_t *k) {
    fmd_kmer_kv_t *retval = NULL;
    if(k) {
        retval = calloc(1, sizeof(fmd_kmer_kv_t));
        retval->str = fmd_string_copy(k->str);
        retval->vid = k->vid;
        retval->offset = k->offset;
    }
    return retval;
}

void fmd_graph_kmer_locations_helper_trav(void *key, void *value, void *params) {
    fmd_vector_t *v = (fmd_vector_t*)params;
    fmd_vector_t *k = value;
    for(int_t i = 0; i < k->size; i++) {
        fmd_kmer_kv_t *kv = k->data[i];
        kv->str = key;
        fmd_vector_append(v, kv);
    }
}

void fmd_graph_kmer_locations(fmd_graph_t *graph, int_t k, fmd_vector_t **kmers, fmd_table_t **kmer_table) {
    int_t V = graph->vertex_list->size;
    fmd_table_t *set;
    fmd_table_init(&set, FMD_HT_INIT_SIZE, &fmd_fstruct_string, &fmd_fstruct_vector);
    for(int_t i = 0; i < V; i++) {
        fmd_vertex_t *vertex = graph->vertex_list->data[i];
        for(int_t j = 0; j < vertex->label->size - k + 1; j++) {
            fmd_string_t *kmer;
            fmd_string_substring(vertex->label, j, j + k, &kmer);
            fmd_kmer_kv_t *loc = calloc(1, sizeof(fmd_kmer_kv_t));
            loc->vid = i;
            loc->offset = j;
            fmd_vector_t *locs;
            void *_;
            bool exists = fmd_table_lookup(set, kmer, &_);
            if(!exists) {
                fmd_vector_init(&locs, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_kmer_kv);
                fmd_vector_append(locs, loc);
                fmd_string_t *set_kmer = fmd_string_copy(kmer);
                fmd_table_insert(set, set_kmer, locs);
            } else {
                locs = (fmd_vector_t*)_;
                fmd_vector_append(locs, loc);
            }
            fmd_string_free(kmer);
        }

        // check the outgoing neighbors
        fmd_vector_t *outgoing;
        fmd_table_lookup(graph->outgoing_neighbors, (void*)vertex->id, &outgoing);
        typedef struct kmer_dfs_ {
            fmd_vertex_t *v;
            fmd_string_t *str;
            int_t remaining;
            vid_t source_vid;
            int_t source_offset;
        } kmer_dfs_t;
        fmd_vector_t *stack;
        fmd_vector_init(&stack, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
        // seed the stack
        for(int_t j = 0; j < outgoing->size; j++) {
            fmd_vertex_t *nvertex = graph->vertex_list->data[(int_t)outgoing->data[j]];
            int_t bound = MIN2(k-1, vertex->label->size);
            for(int_t q = 1; q <= bound; q++) {
                fmd_string_t *suf;
                fmd_string_substring(vertex->label, vertex->label->size-q, vertex->label->size, &suf);
                int_t remaining = k - q;
                kmer_dfs_t *rec = calloc(1, sizeof(kmer_dfs_t));
                rec->remaining = remaining;
                rec->str = suf;
                rec->v = nvertex;
                rec->source_vid = i;
                rec->source_offset = vertex->label->size-q;
                fmd_vector_append(stack, rec);
            }
        }

        while(stack->size) {
            void *_;
            fmd_vector_pop(stack, &_);
            kmer_dfs_t *rec = _;
            int_t n_traversable = MIN2(rec->v->label->size, rec->remaining);
            fmd_string_t *prefix;
            fmd_string_substring(rec->v->label, 0, n_traversable, &prefix);
            fmd_string_concat_mut(rec->str, prefix);
            fmd_string_free(prefix);
            if(n_traversable == rec->remaining) {
                fmd_kmer_kv_t *loc = calloc(1,sizeof(fmd_kmer_kv_t));
                loc->vid = rec->source_vid;
                loc->offset = rec->source_offset;
                bool exists = fmd_table_lookup(set, rec->str, &_);
                fmd_vector_t *locs;
                if(exists) {
                    locs = _;
                    fmd_vector_append(locs,loc);
                } else {
                    fmd_vector_init(&locs, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_kmer_kv);
                    fmd_vector_append(locs, loc);
                    fmd_string_t *set_kmer = fmd_string_copy(rec->str);
                    if(!set_kmer->seq[0]) {
                        printf("wtf\n");
                    }
                    fmd_table_insert(set, set_kmer, locs);
                }
                fmd_string_free(rec->str);
                free(rec);
            } else {
                int_t remaining = rec->remaining - n_traversable;
                fmd_vector_t *outgoingv;
                fmd_table_lookup(graph->outgoing_neighbors, (void*)rec->v->id, &_);
                outgoingv = _;
                if(outgoingv->size) {
                    for (int_t j = 0; j < outgoingv->size; j++) {
                        fmd_vertex_t *nvertex = graph->vertex_list->data[(int_t)outgoingv->data[j]];
                        kmer_dfs_t *nrec = calloc(1, sizeof(kmer_dfs_t));
                        nrec->v = nvertex;
                        nrec->str = fmd_string_copy(rec->str);
                        nrec->remaining = remaining;
                        nrec->source_vid = rec->source_vid;
                        nrec->source_offset = rec->source_offset;
                        fmd_vector_append(stack, nrec);
                    }
                } else { // no outgoing vertices
                    fmd_string_free(rec->str);
                    free(rec);
                }
            }
        }
        fmd_vector_free(stack);
    }
    // gather all strings into a vector and sort
    fmd_vector_t *kmer_vec;
    fmd_vector_init(&kmer_vec, set->size, &fmd_fstruct_kmer_kv);
    fmd_table_traverse(set, kmer_vec, fmd_graph_kmer_locations_helper_trav);
    *kmers = kmer_vec;
    *kmer_table = set;
}
