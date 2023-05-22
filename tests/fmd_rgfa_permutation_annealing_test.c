#include "rgfa_parser.h"
#include "fmd_graph.h"
#include "fmd_constraint_set.h"
#include "fmd_annealing.h"

int main() {
    FILE *rgfa_file = fopen("GRCh38-20-0.10b.gfa", "r");
    rgfa_t *rgfa = rgfa_parse(rgfa_file);
    fmd_graph_t *graph = rgfa_to_fmd_graph(rgfa);

    fmd_vector_t *constraint_sets;
    fmd_constraint_set_enumerate(&constraint_sets, graph, 8);

    fmd_annealing_t *ann;
    fmd_annealing_configure(&ann, graph, constraint_sets, 1000000, 1.5, 0.9999, 10);

    while(fmd_annealing_has_more(ann)) {
        fmd_annealing_iterate(ann);
        if((ann->cur_iter + 1) % 10000 == 0) {
            printf("Iteration %d, best_cost = %lf, temperature = %lf\n", ann->cur_iter + 1, ann->best_cost_so_far, ann->temperature);
        }
    }
    printf("Optimal permutation: \n");
    for(int_t i = 0; i < ann->no_vertices; i++) {
        printf("%ld ", ann->best_permutation_so_far[i]);
    }
    printf("\n");
}

