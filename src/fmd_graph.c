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
}

void fmd_graph_free(fmd_graph_t *graph) {
    fmd_vector_free(graph->vertex_list);
    fmd_table_free(graph->vertices);
    fmd_table_free(graph->incoming_neighbors);
    free(graph);
}

void fmd_graph_insert_vertex(fmd_graph_t *graph, vid_t id, fmd_string_t *label) {
    fmd_vertex_t *vertex;
    fmd_vertex_init(&vertex, id, label);
    fmd_table_insert(graph->vertices, (void*)id, vertex);
    fmd_vector_append(graph->vertex_list, vertex);

    fmd_vector_t *neighbors;
    fmd_vector_init(&neighbors, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_insert(graph->incoming_neighbors, (void*)id, neighbors);
}

void fmd_graph_insert_edge(fmd_graph_t *graph, vid_t source, vid_t destination) {
    fmd_vector_t *neighbors;
    fmd_table_lookup(graph->incoming_neighbors, (void*)destination, &neighbors);
    fmd_vector_append(neighbors, (void*)source);
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
    return copy;
}