#include "gin_gin.h"
#include "gin_constraint_set.h"
#include "gin_annealing.h"

gin_graph_t *graph_1(){
    char *str = "AAGGACTAAGGTAACAAGTAA";
    int_t len = strlen(str);
    gin_string_t *labels[len];

    for (int_t i = 0; i < len; i++) {
        char lab[2];
        lab[0] = str[i];
        lab[1] = 0;
        gin_string_init_cstr(&labels[i], lab);
    }

    gin_graph_t *g;
    gin_graph_init(&g);

    for(int_t i = 0; i < len; i++) {
        gin_graph_insert_vertex(g, i, labels[i]);
    }

    for(int_t i = 0; i < len-1; i++) {
        gin_graph_insert_edge(g, i, i+1);
    }

    return g;
}

int main() {
    gin_graph_t *g = graph_1();

    gin_vector_t *constraints;
    gin_constraint_set_enumerate(&constraints, g, 1, true);

    for(int_t i = 0; i < constraints->size; i++) {
        gin_constraint_set_t *s = constraints->data[i];
        printf("%s : [ ", s->str->seq);
        for(int_t j = 0; j < s->vertices->size; j++) {
            printf("%d ", s->vertices->data[j]);
        }
        printf("]\n");
    }

    gin_annealing_t *ann;
    gin_annealing_configure(&ann, g, NULL, constraints, 1000, 1, 0.99, 1 );

    gin_annealing_iterate_seconds(ann,5);

    srand(42);
    printf("Permutation found with cost %1.lf: [ ", ann->best_cost_so_far);
    for(int_t i = 0; i < ann->no_vertices; i++) {
        printf("%d ", ann->best_permutation_so_far[i]);
    }
    printf("]\n");

    gin_vector_t *perm;
    gin_annealing_get_permutation(ann, &perm);


    // now construct the gin
    printf("Regular:\n");
    gin_gin_t *gin;
    gin_gin_init(&gin, g, NULL, GIN_GIN_DEFAULT_c_0, GIN_GIN_DEFAULT_c_1, 256, 256);


    printf("Optimized:\n");
    gin_gin_t *gin_optim;
    gin_gin_init(&gin_optim, g, perm, GIN_GIN_DEFAULT_c_0, GIN_GIN_DEFAULT_c_1, 256, 256);

    gin_vector_t *paths,*dead;
    int_t forks;
    gin_string_t *query;
    gin_string_init_cstr(&query, "GG");
    gin_gin_stats_t s;
    gin_gin_query_find(gin_optim, NULL, query, -1, &paths, &dead, &s);

    int_t count = 0;
    for(int_t i = 0; i < paths->size; i++) {
        gin_fork_node_t *fork = paths->data[i];
        count += fork->sa_hi - fork->sa_lo;
    }

    printf("no forks optim: %d, count = %d\n", forks, count);

    gin_gin_query_find(gin_optim, NULL, query, -1, &paths, &dead, &s);

    count = 0;
    for(int_t i = 0; i < paths->size; i++) {
        gin_fork_node_t *fork = paths->data[i];
        count += fork->sa_hi - fork->sa_lo;
    }

    printf("no forks regular: %d, count = %d\n", forks, count);


}