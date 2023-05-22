#ifndef FMD_FMD_ANNEALING_H
#define FMD_FMD_ANNEALING_H

#include "fmd_common.h"
#include "fmd_graph.h"

typedef struct fmd_annealing_ {
    byte_t **bin_matrix;
    int_t no_constraints;
    int_t no_vertices;
    vid_t *permutation;
    int_t *block_counts;

    double temperature;
    double cooling_factor;
    double scaling_factor;
    double min_temperature;
    int cur_iter;
} fmd_annealing_t;

void fmd_annealing_configure(fmd_annealing_t **cfg,
                             fmd_graph_t *fmd_graph,
                             fmd_vector_t *constraint_sets,
                             double temperature,
                             double scaling_factor,
                             double cooling_factor,
                             double min_temperature);
void fmd_annealing_launch(fmd_annealing_t *cfg);


#endif //FMD_FMD_ANNEALING_H
