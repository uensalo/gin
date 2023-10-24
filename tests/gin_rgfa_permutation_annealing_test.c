#include "rgfa_parser.h"
#include "gin_graph.h"
#include "gin_constraint_set.h"
#include "gin_annealing.h"
#include "time.h"

double permutation_cost(vid_t *permutation, int_t V, gin_vector_t *constraints) {
    byte_t **M = calloc(V, sizeof(byte_t*)); // store column-wise
    byte_t **Mp = calloc(V, sizeof(byte_t*));
    for(int_t i = 0; i < V; i++){
        M[i] = calloc(constraints->size, sizeof(byte_t));
    }
    for(int_t i = 0; i < constraints->size; i++) {
        gin_constraint_set_t *constraint = constraints->data[i];
        gin_vector_t *constraint_vertices = constraint->vertices;
        for(int_t j = 0; j < constraint_vertices->size; j++) {
            vid_t vertex_no = (vid_t)constraint_vertices->data[j];
            M[vertex_no][i] = 1;
        }
    }
    // permute
    for(int_t i = 0; i < V; i++) {
        Mp[i] = M[(int_t)permutation[i]];
    }
    // compute cost
    double cost = 0;
    for(int_t i = 0; i < constraints->size; i++) {
        bool in_block = false;
        for(int_t j = 0; j < V; j++) {
            if (Mp[j][i] == 1) {
                if (!in_block) {
                    ++cost;
                    in_block = true;
                }
            } else {
                in_block = false;
            }
        }
    }
    return cost;
}

void timestamp() {
    time_t ltime;
    ltime=time(NULL);
    char *str = asctime( localtime(&ltime) );
    str[strlen(str)-1]=0;
    printf("[%s]: ", str);
}

int main() {
    srand(908712892);
    int_t depth = 6;

    FILE *rgfa_file = fopen("GRCh38-20-0.10b.gfa", "r");
    timestamp(); printf("Parsing the graph genome file...\n");
    rgfa_t *rgfa = rgfa_parse(rgfa_file);
    gin_graph_t *graph = rgfa_to_gin_graph(rgfa);
    rgfa_free(rgfa);



    timestamp(); printf("Enumerating constraint sets for depth=%ld\n", depth);
    gin_vector_t *constraint_sets;
    gin_constraint_set_enumerate(&constraint_sets, graph, 4, true);
    timestamp(); printf("Number of constraint sets enumerated: %ld\n", constraint_sets->size);

    timestamp(); printf("Configuring annealing...\n");
    gin_annealing_t *ann;
    gin_annealing_configure(&ann, graph, constraint_sets, NULL, 0, 1.5, 1, 10);
    timestamp(); printf("Annealing configured, running optimization...\n");

    struct timespec t1,t2;
    double time_elapsed = 0;
    double max_time = 1800;
    while(time_elapsed < max_time) {
        clock_gettime(CLOCK_REALTIME, &t1);
        gin_annealing_iterate(ann);
        if((ann->cur_iter + 1) % 10000 == 0) {
            timestamp(); printf("Iteration %d, best_cost = %lf, cur_cost = %lf, temperature = %lf\n", ann->cur_iter + 1, ann->best_cost_so_far, ann->cur_cost, ann->temperature);
        }
        clock_gettime(CLOCK_REALTIME, &t2);
        time_elapsed += (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9;
    }
    printf("Optimal permutation: \n");
    for(int_t i = 0; i < ann->no_vertices; i++) {
        printf("%ld ", ann->best_permutation_so_far[i]);
    }
    printf("\n");
    double pcost = permutation_cost(ann->best_permutation_so_far, ann->no_vertices, constraint_sets);
    printf("Cost comparison: opt = %lf, naive = %lf", ann->best_cost_so_far, pcost);
}

