#ifndef FMD_FMD_GRAPH_H
#define FMD_FMD_GRAPH_H

#include <fmd_string.h>
#include <fmd_table.h>
#include <fmd_tree.h>
#include <fmd_vector.h>

typedef int_t vid_t;
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

typedef struct fmd_graph_flatten_alphabet_ {
    fmd_vector_t *alphabet;
    fmd_table_t *idx2char;
    fmd_table_t *char2idx;
} fmd_graph_flatten_alphabet_t;
void fmd_graph_flatten_alphabet(void *key, void *value, void *params);
typedef struct fmd_graph_ecs_ {
    int_t pos;
    vid_t head_vid;
    vid_t end_vid;
} fmd_graph_ecs_t;
fmd_table_t *fmd_graph_extract_constraint_sets(fmd_graph_t *graph, int_t max_depth);
void fmd_graph_extract_constraint_sets_helper(fmd_vector_t *paths,
                                              fmd_string_t *prefix,
                                              fmd_table_t *constraint_sets,
                                              fmd_vector_t *alphabet,
                                              fmd_table_t *char2idx,
                                              fmd_table_t *idx2char,
                                              fmd_graph_t *graph,
                                              int_t max_depth);

void fmd_graph_flatten_constraints(void *key, void *value, void *p);

#endif //FMD_FMD_GRAPH_H
