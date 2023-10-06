#include "fmd_fmd.h"
#include "fmd_graph.h"

fmd_graph_t *test_graph1() {
    char *l0 = "ACCGTA";
    char *l1 = "ACGTTA";
    char *l2 = "GTTATA";
    char *l3 = "CCGTTA";
    fmd_string_t *l0s, *l1s, *l2s, *l3s;
    fmd_string_init_cstr(&l0s, l0);
    fmd_string_init_cstr(&l1s, l1);
    fmd_string_init_cstr(&l2s, l2);
    fmd_string_init_cstr(&l3s, l3);

    fmd_graph_t *graph;
    fmd_graph_init(&graph);

    fmd_graph_insert_vertex(graph, 0, l0s);
    fmd_graph_insert_vertex(graph, 1, l1s);
    fmd_graph_insert_vertex(graph, 2, l2s);
    fmd_graph_insert_vertex(graph, 3, l3s);

    fmd_graph_insert_edge(graph, 0, 1);
    fmd_graph_insert_edge(graph, 0, 2);
    fmd_graph_insert_edge(graph, 1, 3);
    fmd_graph_insert_edge(graph, 2, 3);

    return graph;
}

int main() {
    fmd_graph_t *g = test_graph1();
    fmd_fmd_t *fmd;
    fmd_fmd_init(&fmd, g, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 1, 1);

    fmd_string_t *query;
    char *querycstr = "TAC";
    fmd_string_init_cstr(&query, querycstr);
    fmd_vector_t *paths, *dead_ends;
    fmd_fmd_stats_t s;
    fmd_fmd_query_find(fmd, NULL, query, -1, &paths, &dead_ends, &s);

    fmd_vector_free(paths);
    fmd_vector_free(dead_ends);
    fmd_fmd_free(fmd);
    fmd_string_free(query);
    fmd_graph_free(g);

    return 0;
}