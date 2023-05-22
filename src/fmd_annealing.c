#include "fmd_annealing.h"

void fmd_annealing_configure(fmd_annealing_t **cfg,
                             fmd_graph_t *graph,
                             fmd_vector_t *constraint_sets,
                             double temperature,
                             double scaling_factor,
                             double cooling_factor,
                             double min_temperature) {
    if(!graph || !constraint_sets) {
        *cfg = NULL;
        return;
    }
    fmd_annealing_t *config = calloc(1, sizeof(fmd_annealing_t));
    if(!config) {
        *cfg = NULL;
        return;
    }
    config->temperature = temperature;
    config->scaling_factor = scaling_factor;
    config->cooling_factor = cooling_factor;
    config->min_temperature = min_temperature;
    config->cur_iter = 0;
    // convert constraint sets into matrices
    config->no_constraints = constraint_sets->size;
    config->no_vertices = graph->vertices->size;

    config->bin_matrix = calloc(config->no_vertices, sizeof(byte_t*)); // store column-wise
    for(int_t i = 0; config->no_constraints; i++){
        config->bin_matrix[i] = calloc(config->no_constraints, sizeof(byte_t));
    }

    config->permutation = calloc(config->no_vertices, sizeof(vid_t));
    for(int_t i = 0; i < config->no_vertices; i++) {
        config->permutation[i] = i; // initial permutation is the identity map
    }

    // now populate the initial constraints
    for(int_t i = 0; i < config->no_constraints; i++) {
        fmd_vector_t *constraint_vertices = constraint_sets->data[i];
        for(int_t j = 0; j < graph->vertices->size; j++) {
            vid_t vertex_no = (vid_t)constraint_vertices->data[j];
            config->bin_matrix[vertex_no][i] = 1;
        }
    }

    // compute and cache run-lengths
    config->block_counts = calloc(config->no_constraints, sizeof(int_t));
    for(int_t i = 0; i < config->no_constraints; i++) {
        int_t cur_run_length = 0;
        for(int_t j = 0; j < config->no_vertices; j++) {

        }
    }
}

void fmd_annealing_launch(fmd_annealing_t *cfg) {

}