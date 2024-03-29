/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_annealing.c is part of gin
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
#include "gin_annealing.h"
#include "math.h"
#include "time.h"

void gin_annealing_configure(gin_annealing_t **cfg,
                             gin_graph_t *graph,
                             gin_vector_t *constraint_sets,
                             gin_vector_t *initial_permutation,
                             double temperature,
                             double scaling_factor,
                             double cooling_factor,
                             double min_temperature) {
    if(!graph || !constraint_sets) {
        *cfg = NULL;
        return;
    }
    gin_annealing_t *config = calloc(1, sizeof(gin_annealing_t));
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
    for(int_t i = 0; i < config->no_vertices; i++){
        config->bin_matrix[i] = calloc(config->no_constraints, sizeof(byte_t));
    }

    config->permutation = calloc(config->no_vertices, sizeof(vid_t));
    vid_t *inverse_permutation = calloc(config->no_vertices, sizeof(vid_t));
    for(int_t i = 0; i < config->no_vertices; i++) {
        config->permutation[i] = initial_permutation ? (vid_t)initial_permutation->data[i] : i;
        inverse_permutation[initial_permutation ? (vid_t)initial_permutation->data[i] : i] = i;
    }
    config->best_permutation_so_far = calloc(config->no_vertices, sizeof(vid_t));
    memcpy(config->best_permutation_so_far, config->permutation, sizeof(vid_t) * config->no_vertices);

    // now populate the initial constraints
#ifdef GIN_OMP
#pragma omp parallel for default(none) shared(config, constraint_sets, inverse_permutation)
#endif
    for(int_t i = 0; i < config->no_constraints; i++) {
        gin_constraint_set_t *constraint = constraint_sets->data[i];
        gin_vector_t *constraint_vertices = constraint->vertices;
        for(int_t j = 0; j < constraint_vertices->size; j++) {
            vid_t permuted_index = inverse_permutation[(vid_t)constraint_vertices->data[j]];
            config->bin_matrix[permuted_index][i] = 1;
        }
    }

    // compute and cache number of runs
    config->cur_cost = 0;
    config->block_counts = calloc(config->no_constraints, sizeof(int_t));
#ifdef GIN_OMP
#pragma omp parallel for default(none) shared(config)
#endif
    for(int_t i = 0; i < config->no_constraints; i++) {
        bool in_block = false;
        for(int_t j = 0; j < config->no_vertices; j++) {
            if (config->bin_matrix[j][i] == 1) {
                if (!in_block) {
#ifndef GIN_OMP
                    ++config->cur_cost;
#endif
                    ++config->block_counts[i];
                    in_block = true;
                }
            } else {
                in_block = false;
            }
        }
    }
#ifdef GIN_OMP
    int_t total_count = 0;
#pragma omp parallel for default(none) shared(config) reduction(+:total_count)
    for(int_t i = 0; i < config->no_constraints; i++) {
        total_count += config->block_counts[i];
    }
    config->cur_cost = (double)total_count;
#endif
    config->best_cost_so_far = config->cur_cost;
    config->next_block_counts = calloc(config->no_constraints, sizeof(int_t));
    *cfg = config;

    free(inverse_permutation);
}

void gin_annealing_step_naive(gin_annealing_t *ann, int_t v1, int_t v2) {
    // swap in the matrix
    byte_t *tmp = ann->bin_matrix[v1];
    ann->bin_matrix[v1] = ann->bin_matrix[v2];
    ann->bin_matrix[v2] = tmp;

    // swap in the permutation
    vid_t tmp_vid = ann->permutation[v1];
    ann->permutation[v1] = ann->permutation[v2];
    ann->permutation[v2] = tmp_vid;

    double next_cost = 0;
    memset(ann->next_block_counts, 0, sizeof(vid_t)*ann->no_constraints);
    for(int_t i = 0; i < ann->no_constraints; i++) {
        bool in_block = false;
        for(int_t j = 0; j < ann->no_vertices; j++) {
            if (ann->bin_matrix[j][i] == 1) {
                if (!in_block) {
                    ++ann->next_block_counts[i];
                    ++next_cost;
                    in_block = true;
                }
            } else {
                in_block = false;
            }
        }
    }
    ann->next_cost = next_cost;
}
void gin_annealing_step_unrolled(gin_annealing_t *ann, int_t v1, int_t v2) {
    ann->next_cost = ann->cur_cost;
    // the only effected runs are at the indices of the corresponding swaps...
    // no need to iterate the whole blocks
#ifdef GIN_OMP
#pragma omp parallel for default(none) shared(ann,v1,v2)
#endif
    for(int_t i = 0; i < ann->no_constraints; i++) {
        if (ann->bin_matrix[v1][i] == ann->bin_matrix[v2][i]) { // no cost change, ignore
            ann->next_block_counts[i] = ann->block_counts[i];
            continue;
        }
        if(v1 - v2 == 1 || v2 - v1 == 1) {
            // special case of adjacent swap
            int_t u = (v1 < v2 ? v1 : v2) - 1;
            int_t d = (v1 < v2 ? v2 : v1) + 1;
            if(u == -1) {
                // adjacent swap on the first row
                if(ann->bin_matrix[u+1][i] == 0) { // ann->bin_matrix[d-1][i] == 1 implied
                    if(ann->bin_matrix[d][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // ann->bin_matrix[d-1][i] == 0 implied
                    if(ann->bin_matrix[d][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else if (d == ann->no_vertices) {
                // adjacent swap on the last row
                if(ann->bin_matrix[u+1][i] == 0) { // ann->bin_matrix[d-1][i] == 1 implied
                    if(ann->bin_matrix[u][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // ann->bin_matrix[d-1][i] == 0 implied
                    if(ann->bin_matrix[u][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else {
                // general case
                if(ann->bin_matrix[u+1][i]==0) { // ann->bin_matrix[d-1][i] == 0 implied
                    if(ann->bin_matrix[u][i] == 1 && ann->bin_matrix[d][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else if(ann->bin_matrix[u][i] == 0 && ann->bin_matrix[d][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // ann->bin_matrix[d-1][i] == 1 implied
                    if(ann->bin_matrix[u][i] == 1 && ann->bin_matrix[d][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else if(ann->bin_matrix[u][i] == 0 && ann->bin_matrix[d][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            }
        } else {
            // case 1: the swap changes the bit
            // handle v1
            if (v1 == 0) {
                // edge case 1: swap is on the first row
                if (ann->bin_matrix[v1][i] == 0) { // implied ann->bin_matrix[v2][i] == 1
                    if (ann->bin_matrix[v1+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v2][i] == 0
                    if (ann->bin_matrix[v1+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else if (v1 == ann->no_vertices-1) {
                // edge case 2: swap is on the last row
                if (ann->bin_matrix[v1][i] == 0) { // implied ann->bin_matrix[v2][i] == 1
                    if (ann->bin_matrix[v1-1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v2][i] == 0
                    if (ann->bin_matrix[v1-1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else {
                // general case
                if (ann->bin_matrix[v1][i] == 0) { // implied ann->bin_matrix[v2][i] == 1
                    if (ann->bin_matrix[v1-1][i] == 1 && ann->bin_matrix[v1+1][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else if (ann->bin_matrix[v1-1][i] == 0 && ann->bin_matrix[v1+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v2][i] == 0
                    if (ann->bin_matrix[v1-1][i] == 1 && ann->bin_matrix[v1+1][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else if (ann->bin_matrix[v1-1][i] == 0 && ann->bin_matrix[v1+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            }
            // handle v2
            if (v2 == 0) {
                // edge case 1: swap is on the first row
                if (ann->bin_matrix[v2][i] == 0) { // implied ann->bin_matrix[v1][i] == 1
                    if (ann->bin_matrix[v2+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v1][i] == 0
                    if (ann->bin_matrix[v2+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else if (v2 == ann->no_vertices-1) {
                // edge case 2: swap is on the last row
                if (ann->bin_matrix[v2][i] == 0) { // implied ann->bin_matrix[v1][i] == 1
                    if (ann->bin_matrix[v2-1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v1][i] == 0
                    if (ann->bin_matrix[v2-1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            } else {
                // general case
                if (ann->bin_matrix[v2][i] == 0) { // implied ann->bin_matrix[v1][i] == 1
                    if (ann->bin_matrix[v2-1][i] == 1 && ann->bin_matrix[v2+1][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else if (ann->bin_matrix[v2-1][i] == 0 && ann->bin_matrix[v2+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                } else { // implied ann->bin_matrix[v1][i] == 0
                    if (ann->bin_matrix[v2-1][i] == 1 && ann->bin_matrix[v2 + 1][i] == 1) {
                        ann->next_block_counts[i] = ann->block_counts[i] + 1;
                        ++ann->next_cost;
                    } else if (ann->bin_matrix[v2-1][i] == 0 && ann->bin_matrix[v2+1][i] == 0) {
                        ann->next_block_counts[i] = ann->block_counts[i] - 1;
                        --ann->next_cost;
                    } else {
                        ann->next_block_counts[i] = ann->block_counts[i];
                    }
                }
            }
        }
    }

    // swap in the matrix
    byte_t *tmp = ann->bin_matrix[v1];
    ann->bin_matrix[v1] = ann->bin_matrix[v2];
    ann->bin_matrix[v2] = tmp;

    // swap in the permutation
    vid_t tmp_vid = ann->permutation[v1];
    ann->permutation[v1] = ann->permutation[v2];
    ann->permutation[v2] = tmp_vid;
}

void gin_annealing_step(gin_annealing_t *ann, int_t v1, int_t v2) {
    byte_t **M = ann->bin_matrix;
    int_t V = ann->no_vertices;
    int_t C = ann->no_constraints;
#ifdef GIN_OMP
#pragma omp parallel for default(none) shared(M,V,C,v1,v2,ann)
#else
    double next_cost = ann->cur_cost;
#endif
    for(int_t i = 0; i < C; i++) {
        if(M[v1][i] == M[v2][i]) { // swap doesn't change anything, move on
            ann->next_block_counts[i] = ann->block_counts[i];
            continue;
        }
        int_t s = v1 < v2 ? v1 : v2;
        int_t b = v1 < v2 ? v2 : v1;
        bool ap = b != s+1;
        byte_t vs = M[s][i]; // vs == 0 implies vb == 1, no need to look the other one up
        int_t a0 = s > 0 ? M[s-1][i] : 0;
        int_t a1 = ap ? M[s+1][i] : 0;
        int_t a2 = ap ? M[b-1][i] : 0;
        int_t a3 = b < V-1 ? M[b+1][i] : 0;
        int_t del = (vs?-1:1)*((a3+a2)-(a1+a0)); // black magic super concise encoding
        ann->next_block_counts[i] = ann->block_counts[i] + del;
#ifndef GIN_OMP
        next_cost += (double)del;
#endif
    }
#ifdef GIN_OMP
    int_t next_count = 0;
#pragma omp parallel for default(none) reduction(+:next_count) shared(ann,C)
    for(int_t i = 0; i < C; i++) {
        next_count += ann->next_block_counts[i];
    }
    double next_cost = (double)next_count;
#endif
    ann->next_cost = next_cost;
    // swap in the matrix
    byte_t *tmp = ann->bin_matrix[v1];
    ann->bin_matrix[v1] = ann->bin_matrix[v2];
    ann->bin_matrix[v2] = tmp;

    // swap in the permutation
    vid_t tmp_vid = ann->permutation[v1];
    ann->permutation[v1] = ann->permutation[v2];
    ann->permutation[v2] = tmp_vid;
}

void gin_annealing_accept(gin_annealing_t *ann) {
    ann->cur_cost = ann->next_cost;
    int_t *tmp = ann->block_counts;
    ann->block_counts = ann->next_block_counts;
    ann->next_block_counts = tmp;
}

void gin_annealing_reject(gin_annealing_t *ann, int_t v1, int_t v2) {
    // swap back in the matrix
    byte_t *tmp = ann->bin_matrix[v1];
    ann->bin_matrix[v1] = ann->bin_matrix[v2];
    ann->bin_matrix[v2] = tmp;

    // swap back in the permutation
    vid_t tmp_vid = ann->permutation[v1];
    ann->permutation[v1] = ann->permutation[v2];
    ann->permutation[v2] = tmp_vid;
}

bool gin_annealing_has_more(gin_annealing_t *ann) {
    return ann->temperature >= ann->min_temperature;
}

void gin_annealing_iterate(gin_annealing_t *ann) {
    int_t v1 = rand() % ann->no_vertices;
    int_t v2;
    do {
        v2 = rand() % ann->no_vertices;
    } while (v1 == v2);

    gin_annealing_step(ann, v1, v2);
    double acceptance_prob = ann->next_cost < ann->cur_cost ? 1.0 :
                             exp((ann->cur_cost - ann->next_cost) / (ann->temperature * ann->scaling_factor));

    // sample
    if (acceptance_prob < ((double)rand() / (double)RAND_MAX)) {
        gin_annealing_reject(ann, v1, v2);
    } else {
        gin_annealing_accept(ann);
    }

    if(ann->cur_cost < ann->best_cost_so_far) {
        ann->best_cost_so_far = ann->cur_cost;
        memcpy(ann->best_permutation_so_far, ann->permutation, sizeof(vid_t) * ann->no_vertices);
    }

    ann->temperature *= ann->cooling_factor;
    ann->cur_iter += 1;
}

void gin_annealing_iterate_until_end(gin_annealing_t *ann) {
    while(ann->temperature >= ann->min_temperature) {
        int_t v1 = rand() % ann->no_vertices;
        int_t v2;
        do {
            v2 = rand() % ann->no_vertices;
        } while (v1 == v2);

        gin_annealing_step(ann, v1, v2);

        double acceptance_prob = ann->next_cost < ann->cur_cost ? 1.0 :
                                 exp((ann->cur_cost - ann->next_cost) / (ann->temperature * ann->scaling_factor));

        // sample
        if (acceptance_prob < ((double)rand() / (double)RAND_MAX)) {
            gin_annealing_reject(ann, v1, v2);
        } else {
            gin_annealing_accept(ann);
        }
        ann->temperature *= ann->cooling_factor;
        ++ann->cur_iter;
    }
}

void gin_annealing_iterate_seconds(gin_annealing_t *ann, int_t seconds) {
    struct timespec t1,t2;
    double time_elapsed = 0;
    double max_time = (double)seconds;
    while(time_elapsed < max_time) {
        clock_gettime(CLOCK_REALTIME, &t1);
        gin_annealing_iterate(ann);
        clock_gettime(CLOCK_REALTIME, &t2);
        time_elapsed += (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9;
    }
}

void gin_annealing_iterate_seconds_verbose(gin_annealing_t *ann, int_t seconds) {
    struct timespec t1,t2;
    double time_elapsed = 0;
    double max_time = (double)seconds;
    while(time_elapsed < max_time) {
        clock_gettime(CLOCK_REALTIME, &t1);
        gin_annealing_iterate(ann);
        clock_gettime(CLOCK_REALTIME, &t2);
        if((ann->cur_iter + 1) % 10000000 == 0) {
            printf("Iteration %d, best_cost = %lf, cur_cost = %lf, temperature = %lf\n", ann->cur_iter + 1, ann->best_cost_so_far, ann->cur_cost, ann->temperature);
        }
        time_elapsed += (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9;
    }
}

void gin_annealing_get_permutation(gin_annealing_t *ann, gin_vector_t **permutation) {
    gin_vector_t *perm;
    gin_vector_init(&perm, ann->no_vertices, &prm_fstruct);
    for(int_t i = 0; i < ann->no_vertices; i++) {
        gin_vector_append(perm, (void*)ann->best_permutation_so_far[i]);
    }
    *permutation = perm;
}
