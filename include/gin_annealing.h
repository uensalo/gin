/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_annealing.h is part of gin
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
#ifndef GIN_GIN_ANNEALING_H
#define GIN_GIN_ANNEALING_H

#include "gin_common.h"
#include "gin_graph.h"
#include "gin_constraint_set.h"

typedef struct gin_annealing_ {
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
} gin_annealing_t;

void gin_annealing_step_naive(gin_annealing_t *ann, int_t v1, int_t v2);
void gin_annealing_step_unrolled(gin_annealing_t *ann, int_t v1, int_t v2);
void gin_annealing_step(gin_annealing_t *ann, int_t v1, int_t v2);
void gin_annealing_accept(gin_annealing_t *ann);
void gin_annealing_reject(gin_annealing_t *ann, int_t v1, int_t v2);

void gin_annealing_configure(gin_annealing_t **cfg,
                             gin_graph_t *gin_graph,
                             gin_vector_t *constraint_sets,
                             gin_vector_t *initial_permutation,
                             double temperature,
                             double scaling_factor,
                             double cooling_factor,
                             double min_temperature);


bool gin_annealing_has_more(gin_annealing_t *ann);
void gin_annealing_iterate(gin_annealing_t *ann);
void gin_annealing_iterate_until_end(gin_annealing_t *ann);
void gin_annealing_iterate_seconds(gin_annealing_t *ann, int_t seconds);
void gin_annealing_iterate_seconds_verbose(gin_annealing_t *ann, int_t seconds);

void gin_annealing_get_permutation(gin_annealing_t *ann, gin_vector_t **permutation);



#endif //GIN_GIN_ANNEALING_H
