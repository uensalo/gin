#include "fmd_fmd.h"
#include "fmd_constraint_set.h"
#include "fmd_annealing.h"

fmd_graph_t *graph_1(){
    char *str = "AAGGACTAAGGTAACAAGTAA";
    int_t len = strlen(str);
    fmd_string_t *labels[len];

    for (int_t i = 0; i < len; i++) {
        char lab[2];
        lab[0] = str[i];
        lab[1] = 0;
        fmd_string_init_cstr(&labels[i], lab);
    }

    fmd_graph_t *g;
    fmd_graph_init(&g);

    for(int_t i = 0; i < len; i++) {
        fmd_graph_insert_vertex(g, i, labels[i]);
    }

    for(int_t i = 0; i < len-1; i++) {
        fmd_graph_insert_edge(g, i, i+1);
    }

    return g;
}

int main() {
    fmd_graph_t *g = graph_1();

    fmd_vector_t *constraints;
    fmd_constraint_set_enumerate(&constraints, g, 1, true);

    for(int_t i = 0; i < constraints->size; i++) {
        fmd_constraint_set_t *s = constraints->data[i];
        printf("%s : [ ", s->str->seq);
        for(int_t j = 0; j < s->vertices->size; j++) {
            printf("%d ", s->vertices->data[j]);
        }
        printf("]\n");
    }

    fmd_annealing_t *ann;
    fmd_annealing_configure(&ann, g, NULL, constraints, 1000, 1, 0.99, 1 );

    fmd_annealing_iterate_seconds(ann,5);

    srand(42);
    printf("Permutation found with cost %1.lf: [ ", ann->best_cost_so_far);
    for(int_t i = 0; i < ann->no_vertices; i++) {
        printf("%d ", ann->best_permutation_so_far[i]);
    }
    printf("]\n");

    fmd_vector_t *perm;
    fmd_annealing_get_permutation(ann, &perm);


    // now construct the fmd
    printf("Regular:\n");
    fmd_fmd_t *fmd;
    fmd_fmd_init(&fmd, g, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 256, 256);


    printf("Optimized:\n");
    fmd_fmd_t *fmd_optim;
    fmd_fmd_init(&fmd_optim, g, perm, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 256, 256);

    fmd_vector_t *paths,*dead;
    int_t forks;
    fmd_string_t *query;
    fmd_string_init_cstr(&query, "GG");
    fmd_fmd_stats_t s;
    fmd_fmd_query_find(fmd_optim, NULL, query, -1, &paths, &dead, &s);

    int_t count = 0;
    for(int_t i = 0; i < paths->size; i++) {
        fmd_fork_node_t *fork = paths->data[i];
        count += fork->sa_hi - fork->sa_lo;
    }

    printf("no forks optim: %d, count = %d\n", forks, count);

    fmd_fmd_query_find(fmd_optim, NULL, query, -1, &paths, &dead, &s);

    count = 0;
    for(int_t i = 0; i < paths->size; i++) {
        fmd_fork_node_t *fork = paths->data[i];
        count += fork->sa_hi - fork->sa_lo;
    }

    printf("no forks regular: %d, count = %d\n", forks, count);


}