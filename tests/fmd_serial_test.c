#include "fmd_graph.h"
#include "fmd_string.h"
#include "fmd_vector.h"
#include "fmd_fmd.h"

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
    printf("================================================================================\n");
    printf("= GRAPH INITIALIZATION TEST START                                              =\n");
    printf("================================================================================\n");
    fmd_graph_t *g1 = test_graph1();
    for(int_t i = 0; i < 4; i++) {
        fmd_vector_t *incoming;
        fmd_table_lookup(g1->incoming_neighbors, (void*)i, (void*)&incoming);
        printf("Vertex %lld has label: %s", i, ((fmd_vertex_t*)g1->vertex_list->data[i])->label->seq);
        if(incoming) {
            printf(" and incoming neighbors: { ");
            for (int_t k = 0; k < incoming->size; k++) {
                printf("%d ", incoming->data[k]);
            }
            printf("}\n");
        } else {
            printf(" and incoming neighbors: { }\n");
        }
    }
    printf("================================================================================\n");
    printf("= GRAPH INITIALIZATION TEST END                                                =\n");
    printf("================================================================================\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST START                                       =\n");
    printf("= No permutations, fixed binary permutation encoding                           =\n");
    printf("================================================================================\n");
    fmd_fmd_t *fmd;
    fmd_fmd_init(&fmd, g1, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 1, 1);
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST END                                         =\n");
    printf("================================================================================\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST START                                         =\n");
    printf("= No permutations, fixed binary permutation encoding, serial                   =\n");
    printf("================================================================================\n");
    fmd_string_t *query1;
    char *query1cstr = "GTTA";
    fmd_string_init_cstr(&query1, query1cstr);
    //count_t count = fmd_fmd_query_count(fmd, query1);
    fmd_vector_t *locs = fmd_fmd_query_locate(fmd, query1);

    //printf("Count for %s : %lld\n", query1cstr, count);
    printf("Locate for %s : { ", query1cstr);
    for(int_t i = 0; i < locs->size; i++) {
        printf("%d ", locs->data[i]);
    }
    printf("}\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST END                                           =\n");
    printf("================================================================================\n");
}