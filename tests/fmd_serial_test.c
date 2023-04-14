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

fmd_graph_t *test_graph2() {
    char *l0 = "AACG";
    char *l1 = "GGTA";
    char *l2 = "CGAA";
    char *l3 = "TTGATT";
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
    fmd_graph_insert_edge(graph, 1, 2);
    fmd_graph_insert_edge(graph, 2, 0);
    fmd_graph_insert_edge(graph, 2, 3);

    return graph;
}

fmd_graph_t *test_graph3() {
    char *l0 = "AAAA";
    char *l1 = "CCCC";
    char *l2 = "GGGG";
    char *l3 = "TTTT";
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
    fmd_graph_insert_edge(graph, 1, 0);

    fmd_graph_insert_edge(graph, 2, 3);
    fmd_graph_insert_edge(graph, 3, 2);

    fmd_graph_insert_edge(graph, 1, 2);
    fmd_graph_insert_edge(graph, 2, 1);

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
    char *query1cstr = "TAC";
    fmd_string_init_cstr(&query1, query1cstr);
    //count_t count = fmd_fmd_query_count(fmd, query1);
    //fmd_vector_t *locs = fmd_fmd_query_locate_basic(fmd, query1);
    //printf("Locate for %s : { ", query1cstr);
    //printf("}\n");
    fmd_vector_t *paths,*dead;
    fmd_fmd_query_locate_paths(fmd, query1, &paths, &dead);

    //printf("Count for %s : %lld\n", query1cstr, count);
    int_t total_count = 0;
    printf("Valid paths for: %s\n", query1cstr);
    for(int_t i = 0; i < paths->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)paths->data[i];
        total_count += cur->sa_hi - cur->sa_lo;
        printf("Alive fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Dead explorations and partial matches: \n");
    for(int_t i = 0; i < dead->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)dead->data[i];
        printf("Dead fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Actual count: %d\n", total_count);

    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST END                                           =\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST START                                       =\n");
    printf("= Cyclic graph, no permutations, fixed binary permutation encoding             =\n");
    printf("================================================================================\n");
    fmd_fmd_t *fmd_cyclic;
    fmd_graph_t *g2 = test_graph2();
    fmd_fmd_init(&fmd_cyclic, g2, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 1, 1);
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST END                                         =\n");
    printf("================================================================================\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST START                                         =\n");
    printf("= Cyclic graph, no permutations, fixed binary permutation encoding             =\n");
    printf("================================================================================\n");
    fmd_string_t *query2;
    char *query2cstr = "AACGGGTACGAATTGATT";
    fmd_string_init_cstr(&query2, query2cstr);
    fmd_vector_t *paths_cyclic,*dead_cyclic;
    fmd_fmd_query_locate_paths(fmd_cyclic, query2, &paths_cyclic, &dead_cyclic);
    int_t total_count_cycle = 0;
    printf("Valid paths for: %s\n", query2cstr);
    for(int_t i = 0; i < paths_cyclic->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)paths_cyclic->data[i];
        total_count_cycle += cur->sa_hi - cur->sa_lo;
        printf("Alive fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Dead explorations and partial matches: \n");
    for(int_t i = 0; i < dead_cyclic->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)dead_cyclic->data[i];
        printf("Dead fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Actual count: %d\n", total_count_cycle);

    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST END                                           =\n");
    printf("================================================================================\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST START                                       =\n");
    printf("= Compressed DFA, no permutations, fixed binary permutation encoding           =\n");
    printf("================================================================================\n");
    fmd_fmd_t *fmd_cyclic2;
    fmd_graph_t *g3 = test_graph3();
    fmd_fmd_init(&fmd_cyclic2, g3, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, 1, 1);
    printf("================================================================================\n");
    printf("= FM-DIRECTORY INITIALIZATION TEST END                                         =\n");
    printf("================================================================================\n");
    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST START                                         =\n");
    printf("= Compressed DFA, no permutations, fixed binary permutation encoding           =\n");
    printf("================================================================================\n");
    fmd_string_t *query3;
    char *query3cstr = "CCAAAACCCCGGGGTTTTGGGGCCCCA";
    fmd_string_init_cstr(&query3, query3cstr);
    fmd_vector_t *paths_cyclic2,*dead_cyclic2;
    fmd_fmd_query_locate_paths(fmd_cyclic2, query3, &paths_cyclic2, &dead_cyclic2);
    int_t total_count_cycle3 = 0;
    printf("Valid paths for: %s\n", query3cstr);
    for(int_t i = 0; i < paths_cyclic2->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)paths_cyclic2->data[i];
        total_count_cycle3 += cur->sa_hi - cur->sa_lo;
        printf("Alive fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Dead explorations and partial matches: \n");
    for(int_t i = 0; i < dead_cyclic2->size; i++) {
        fmd_fork_node_t *cur = (fmd_fork_node_t*)dead_cyclic2->data[i];
        printf("Dead fork leaf node %d, sa[lo,hi) = [%d,%d)\n", i, cur->sa_lo, cur->sa_hi);
        while(cur) {
            printf("\tpos, [vlo,vhi) = %d, [%d,%d)\n", cur->pos, cur->vertex_lo, cur->vertex_hi);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    printf("Actual count: %d\n", total_count_cycle3);

    printf("================================================================================\n");
    printf("= FM-DIRECTORY COUNT-LOCATE TEST END                                           =\n");
    printf("================================================================================\n");
}