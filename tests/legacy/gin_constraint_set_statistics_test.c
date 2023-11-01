#include "gin_constraint_set.h"
#include <time.h>

void timestamp() {
    time_t ltime;
    ltime=time(NULL);
    char *str = asctime( localtime(&ltime) );
    str[strlen(str)-1]=0;
    printf("[%s]: ", str);
}

gin_graph_t *generate_random_graph(int no_vertices,
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

    gin_graph_t *graph;
    gin_graph_init(&graph);

    for(int i = 0; i < no_vertices; i++) {
        gin_string_t *label;
        gin_string_init_cstr(&label, labels[i]);
        gin_graph_insert_vertex(graph, i , label);
        free(labels[i]);
    }
    free(labels);

    gin_table_t *edges;
    gin_table_init(&edges, GIN_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
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
            if(!gin_table_lookup(edges,(void*)pairing,&_)) {
                gin_table_insert(edges, (void*)pairing, NULL);
                gin_graph_insert_edge(graph, src, dst);
                break;
            }
        }
    }
    return graph;
}

int main() {
    srand(1234567);
    int_t depth = 10;
    gin_graph_t *graph = generate_random_graph(1000, 3000, 5, 6);
    gin_vector_t *constraints;
    gin_constraint_set_enumerate(&constraints, graph, depth);

    int_t cardinality[depth];
    int_t count[depth];
    for(int_t i = 0; i < depth; i++) {
        cardinality[i] = count[i] = 0;
    }
    for(int_t i = 0; i < constraints->size; i++) {
        gin_constraint_set_t *constraint = constraints->data[i];
        printf("%s: [ ", constraint->str->seq);
        for(int j = 0; j < constraint->vertices->size; j++) {
            printf("%ld ", (int_t)constraint->vertices->data[j]);
        }
        printf("]\n");
        cardinality[constraint->str->size-1] += constraint->vertices->size;
        ++count[constraint->str->size-1];
    }
    for(int_t i = 0; i < depth; i++) {
        printf("%d\t%d\t%d\t%4.4lf\n", i+1, cardinality[i], count[i], cardinality[i] / (double)count[i]);
    }

    return 0;
}