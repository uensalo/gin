#ifndef FMD_FMD_ANNEALING_H
#define FMD_FMD_ANNEALING_H

#include "fmd_common.h"
#include "fmd_graph.h"
#include "fmd_constraint_set.h"

typedef struct fmd_annealing_ {
    byte_t **bin_matrix;
    int_t no_constraints;
    int_t no_vertices;
    vid_t *permutation;
    int_t *block_counts;
    int_t *next_block_counts;

    double temperature;
    double cooling_factor;
    double scaling_factor;
    double min_temperature;

    double cur_cost;
    int cur_iter;
    double next_cost;

    vid_t *best_permutation_so_far;
    double best_cost_so_far;
} fmd_annealing_t;

void fmd_annealing_step_naive(fmd_annealing_t *ann, int_t v1, int_t v2);
void fmd_annealing_step_unrolled(fmd_annealing_t *ann, int_t v1, int_t v2);
void fmd_annealing_step(fmd_annealing_t *ann, int_t v1, int_t v2);
void fmd_annealing_accept(fmd_annealing_t *ann);
void fmd_annealing_reject(fmd_annealing_t *ann, int_t v1, int_t v2);

void fmd_annealing_configure(fmd_annealing_t **cfg,
                             fmd_graph_t *fmd_graph,
                             fmd_vector_t *constraint_sets,
                             double temperature,
                             double scaling_factor,
                             double cooling_factor,
                             double min_temperature);


bool fmd_annealing_has_more(fmd_annealing_t *ann);
void fmd_annealing_iterate(fmd_annealing_t *ann);
void fmd_annealing_iterate_until_end(fmd_annealing_t *ann);
void fmd_annealing_iterate_seconds(fmd_annealing_t *ann, int_t seconds);
void fmd_annealing_iterate_seconds_verbose(fmd_annealing_t *ann, int_t seconds);


#endif //FMD_FMD_ANNEALING_H
