#ifndef FMD_FMD_GRAPH_H
#define FMD_FMD_GRAPH_H

#include <fmd_string.h>
#include <fmd_table.h>
#include <fmd_tree.h>
#include <fmd_vector.h>

typedef int vid_t;
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
    fmd_vertex_comp,
    fmd_vertex_hash,
    fmd_vertex_free,
    fmd_vertex_copy
};

typedef struct {
    fmd_vector_t *vertex_list;
    fmd_table_t *vertices;
    fmd_table_t *incoming_neighbors;
} fmd_graph_t;

// Graph functions
void fmd_graph_init(fmd_graph_t **graph);
void fmd_graph_free(fmd_graph_t *graph);
void fmd_graph_insert_vertex(fmd_graph_t *graph, vid_t id, fmd_string_t *label);
void fmd_graph_insert_edge(fmd_graph_t *graph, int source, int destination);
uint_t fmd_graph_hash(fmd_graph_t *graph);
int fmd_graph_comp(fmd_graph_t *g1, fmd_graph_t *g2);
fmd_graph_t *fmd_graph_copy(fmd_graph_t *graph);

static fmd_fstruct_t fmd_fstruct_graph = {
    fmd_graph_comp,
    fmd_graph_hash,
    fmd_graph_free,
    fmd_graph_copy
};

#endif //FMD_FMD_GRAPH_H
