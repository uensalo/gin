#include "gin_gin.h"
#include "gin_graph.h"

gin_graph_t *test_graph1() {
    char *l0 = "ACCGTA";
    char *l1 = "ACGTTA";
    char *l2 = "GTTATA";
    char *l3 = "CCGTTA";
    gin_string_t *l0s, *l1s, *l2s, *l3s;
    gin_string_init_cstr(&l0s, l0);
    gin_string_init_cstr(&l1s, l1);
    gin_string_init_cstr(&l2s, l2);
    gin_string_init_cstr(&l3s, l3);

    gin_graph_t *graph;
    gin_graph_init(&graph);

    gin_graph_insert_vertex(graph, 0, l0s);
    gin_graph_insert_vertex(graph, 1, l1s);
    gin_graph_insert_vertex(graph, 2, l2s);
    gin_graph_insert_vertex(graph, 3, l3s);

    gin_graph_insert_edge(graph, 0, 1);
    gin_graph_insert_edge(graph, 0, 2);
    gin_graph_insert_edge(graph, 1, 3);
    gin_graph_insert_edge(graph, 2, 3);

    return graph;
}

int main() {
    gin_graph_t *g = test_graph1();
    gin_gin_t *gin;
    gin_gin_init(&gin, g, NULL, GIN_GIN_DEFAULT_c_0, GIN_GIN_DEFAULT_c_1, 1, 1);

    gin_string_t *query;
    char *querycstr = "TAC";
    gin_string_init_cstr(&query, querycstr);
    gin_vector_t *paths, *dead_ends;
    gin_gin_stats_t s;
    gin_gin_query_find(gin, NULL, query, -1, &paths, &dead_ends, &s);

    gin_vector_free(paths);
    gin_vector_free(dead_ends);
    gin_gin_free(gin);
    gin_string_free(query);
    gin_graph_free(g);

    return 0;
}