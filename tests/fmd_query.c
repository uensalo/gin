#include "fmd_fmd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void timestamp() {
    time_t ltime;
    ltime=time(NULL);
    char *str = asctime( localtime(&ltime) );
    str[strlen(str)-1]=0;
    printf("[%s]: ", str);
}

fmd_graph_t *generate_random_graph(int no_vertices,
                                   int no_edges,
                                   int min_label_len, int max_label_len) {
    char charset[] = "ACGT";
    int charsetlen = 4;
    // generate uniformly random labels of random length
    char **labels = calloc(no_vertices, sizeof(char*));
    for(int i = 0; i < no_vertices; i++) {
        int label_len = rand() % (max_label_len - min_label_len) + min_label_len;
        labels[i] = calloc(label_len+1, sizeof(char));
        for(int j = 0; j < label_len; j++) {
            labels[i][j] = charset[rand() % charsetlen];
        }
    }

    fmd_graph_t *graph;
    fmd_graph_init(&graph);

    for(int i = 0; i < no_vertices; i++) {
        fmd_string_t *label;
        fmd_string_init_cstr(&label, labels[i]);
        fmd_graph_insert_vertex(graph, i , label);
        free(labels[i]);
    }
    free(labels);

    fmd_table_t *edges;
    fmd_table_init(&edges, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    for(int j = 0; j < no_edges; j++) {
        while(1) {
            int_t src = rand() % no_vertices;
            int_t dst = rand() % no_vertices;
            while (src == dst) {
                dst = rand() % no_vertices;
            }
            int_t *_;
            int_t sum = src + dst;
            int_t pairing =  (sum * (sum + 1)) / 2 + dst;
            if(!fmd_table_lookup(edges,(void*)pairing,&_)) {
                fmd_table_insert(edges, (void*)pairing, NULL);
                fmd_graph_insert_edge(graph, src, dst);
                break;
            }
        }
    }
    return graph;
}

int test(int no_vertices, int no_edges, int minl, int maxl, int rank_sample_rate, int isa_sample_rate, int no_queries, int minq, int maxq, int sup_print) {
    if(!sup_print) {
        double avg_enc_len = (double)((double)no_vertices * (minl + maxl)/2 + no_vertices * 2 + 1 + (double)no_vertices * fmd_ceil_log2(no_vertices))/1000000.0;
        printf("================================================================================\n");
        printf("TEST: V=%d, E=%d, minl=%d, maxl=%d, B=%d, ISA=%d, AVGL=%lf MBP\n",no_vertices,no_edges,minl,maxl,rank_sample_rate,isa_sample_rate, avg_enc_len);
        printf("================================================================================\n");
        timestamp();printf("Generating random string graph\n");
    }
    fmd_graph_t *random_graph = generate_random_graph(no_vertices, no_edges, minl, maxl);
    if(!sup_print) {
        timestamp();printf("Generated. Constructing FMD\n");
    }
    fmd_fmd_t *fmd;
    fmd_fmd_init(&fmd, random_graph, NULL, FMD_FMD_DEFAULT_c_0, FMD_FMD_DEFAULT_c_1, rank_sample_rate, isa_sample_rate);

    if(!sup_print) {
        timestamp();printf("FMD constructed, generating random queries\n");
    }

    char charset[] = "ACGT";
    int charsetlen = 4;
    struct timespec t1,t2;
    // generate uniformly random labels of random length
    double time_elapsed = 0;
    for(int i = 0; i < no_queries; i++) {
        int qlen = rand() % (maxq - minq) + minq;
        char *querycstr = calloc(qlen+1, sizeof(char));
        for(int j = 0; j < qlen; j++) {
            querycstr[j] = charset[rand() % charsetlen];
        }
        fmd_string_t *query;
        fmd_string_init_cstr(&query, querycstr);
        free(querycstr);
        fmd_vector_t *found, *dead_ends;
        clock_gettime(CLOCK_REALTIME, &t1);
        fmd_fmd_query_locate_paths(fmd, query, &found, &dead_ends);
        clock_gettime(CLOCK_REALTIME, &t2);
        time_elapsed += (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) / 1e9;
        fmd_fmd_locate_paths_result_free(found, dead_ends);
        fmd_string_free(query);
    }

    if(!sup_print) {
        timestamp();printf("Time elapsed: %.6lf, Per query: %.6lf...\n", time_elapsed, time_elapsed / no_queries);
        timestamp();printf("Cleaning up...\n");
    }
    fmd_graph_free(random_graph);
    fmd_fmd_free(fmd);
    return 0;
}

int main() {
    srand(143514355);
    int no_random_tests = 1;

    int minV = 10000;
    int maxV = 10001;

    int noQ  = 1000;
    int minq = 10;
    int maxq = 300;

    int minE = 100000;
    int maxE = 100001;

    int min_minl = 10000;
    int max_minl = 200000;

    int min_maxl = 10000;
    int max_maxl = 200000;

    int min_rank = 63;
    int max_rank = 65;

    int min_isa = 63;
    int max_isa = 64;

    timestamp(); printf("Tests begin.\n");
    for (int i = 0; i < no_random_tests; i++) {

        timestamp();
        printf("%d/%d complete.\n", i, no_random_tests);

        int V;
        int E;
        do {
            V = rand() % (maxV - minV) + minV;
            E = rand() % (maxE - minE) + minE;
        } while(E > V * V / 2);
        int minl, maxl;
        do {
            minl = rand() % (max_minl - min_minl) + min_minl;
            maxl = rand() % (max_maxl - min_maxl) + min_maxl;
        } while (minl >= maxl);
        int B = rand() % (max_rank - min_rank) + min_rank;
        int ISA = rand() % (max_isa - min_isa) + min_isa;

        int result = test(V, E, minl, maxl, B, ISA, noQ, minq, maxq, 0);
    }
    return 0;
}