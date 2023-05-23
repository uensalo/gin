#include "fmd_annealing.h"
#include "math.h"

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
    for(int_t i = 0; i < config->no_vertices; i++){
        config->bin_matrix[i] = calloc(config->no_constraints, sizeof(byte_t));
    }

    config->permutation = calloc(config->no_vertices, sizeof(vid_t));
    for(int_t i = 0; i < config->no_vertices; i++) {
        config->permutation[i] = i; // initial permutation is the identity map
    }
    config->best_permutation_so_far = calloc(config->no_vertices, sizeof(vid_t));

    // now populate the initial constraints
#pragma omp parallel for default(none) shared(config, constraint_sets)
    for(int_t i = 0; i < config->no_constraints; i++) {
        fmd_constraint_set_t *constraint = constraint_sets->data[i];
        fmd_vector_t *constraint_vertices = constraint->vertices;
        for(int_t j = 0; j < constraint_vertices->size; j++) {
            vid_t vertex_no = (vid_t)constraint_vertices->data[j];
            config->bin_matrix[vertex_no][i] = 1;
        }
    }

    // compute and cache number of runs
    config->cur_cost = 0;
    config->block_counts = calloc(config->no_constraints, sizeof(int_t));
    for(int_t i = 0; i < config->no_constraints; i++) {
        bool in_block = false;
        for(int_t j = 0; j < config->no_vertices; j++) {
            if (config->bin_matrix[j][i] == 1) {
                if (!in_block) {
                    ++config->cur_cost;
                    ++config->block_counts[i];
                    in_block = true;
                }
            } else {
                in_block = false;
            }
        }
    }
    config->best_cost_so_far = config->cur_cost;
    config->next_block_counts = calloc(config->no_constraints, sizeof(int_t));
    *cfg = config;
}

void fmd_annealing_step_naive(fmd_annealing_t *ann, int_t v1, int_t v2) {
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

void fmd_annealing_step(fmd_annealing_t *ann, int_t v1, int_t v2) {
    ann->next_cost = ann->cur_cost;
    // the only effected runs are at the indices of the corresponding swaps...
    // no need to iterate the whole blocks
#pragma omp parallel for default(none) shared(ann,v1,v2)
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
                    if (ann->bin_matrix[v2-1][i] == 1 && ann->bin_matrix[v2 + 1][i] == 1) {
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

void fmd_annealing_accept(fmd_annealing_t *ann) {
    ann->cur_cost = ann->next_cost;
    ann->block_counts = ann->next_block_counts;
}

void fmd_annealing_reject(fmd_annealing_t *ann, int_t v1, int_t v2) {
    // swap back in the matrix
    byte_t *tmp = ann->bin_matrix[v1];
    ann->bin_matrix[v1] = ann->bin_matrix[v2];
    ann->bin_matrix[v2] = tmp;

    // swap back in the permutation
    vid_t tmp_vid = ann->permutation[v1];
    ann->permutation[v1] = ann->permutation[v2];
    ann->permutation[v2] = tmp_vid;
}

bool fmd_annealing_has_more(fmd_annealing_t *ann) {
    return ann->temperature >= ann->min_temperature;
}

void fmd_annealing_iterate(fmd_annealing_t *ann) {
    int_t v1 = rand() % ann->no_vertices;
    int_t v2;
    do {
        v2 = rand() % ann->no_vertices;
    } while (v1 == v2);

    fmd_annealing_step(ann, v1, v2);
    // DEBUG
    /*
    double opt_cost = ann->next_cost;
    fmd_annealing_step_naive(ann,v1,v2);
    double naive_cost = ann->next_cost;
    if(opt_cost != naive_cost) {
        fprintf(stderr, "costs are not the same with diff opt-naive=%lf\n", opt_cost - naive_cost);
        //exit(-1);
    }
     */
    // DEBUG

    double acceptance_prob = ann->next_cost < ann->cur_cost ? 1.0 :
                             exp((ann->cur_cost - ann->next_cost) / (ann->temperature * ann->scaling_factor));

    // sample
    if (acceptance_prob < ((double)rand() / (double)RAND_MAX)) {
        fmd_annealing_reject(ann, v1, v2);
    } else {
        fmd_annealing_accept(ann);
    }

    if(ann->cur_cost < ann->best_cost_so_far) {
        ann->best_cost_so_far = ann->cur_cost;
        memcpy(ann->best_permutation_so_far, ann->permutation, sizeof(vid_t) * ann->no_vertices);
    }

    ann->temperature *= ann->cooling_factor;
    ann->cur_iter += 1;
}

void fmd_annealing_iterate_until_end(fmd_annealing_t *ann) {
    while(ann->temperature >= ann->min_temperature) {
        int_t v1 = rand() % ann->no_vertices;
        int_t v2;
        do {
            v2 = rand() % ann->no_vertices;
        } while (v1 == v2);

        fmd_annealing_step(ann, v1, v2);

        double acceptance_prob = ann->next_cost < ann->cur_cost ? 1.0 :
                                 exp((ann->cur_cost - ann->next_cost) / (ann->temperature * ann->scaling_factor));

        // sample
        if (acceptance_prob < ((double)rand() / (double)RAND_MAX)) {
            fmd_annealing_reject(ann, v1, v2);
        } else {
            fmd_annealing_accept(ann);
        }
        ann->temperature *= ann->cooling_factor;
        ann->cur_iter += 1;
    }
}