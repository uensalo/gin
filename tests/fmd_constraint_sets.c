#include "fmd_graph.h"
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

typedef struct cs_ {
    fmd_string_t *str;
    fmd_vector_t *vertices;
} cs_t;

void cs_free(cs_t *cs) {
    free(cs);
}
uint_t cs_hash(cs_t *c1) {
    return fmd_string_hash(c1->str);
}
int cs_comp(cs_t *c1, cs_t *c2) {
    if (c1->str->size < c2->str->size) {
        return -1;
    } else if (c1->str->size > c2->str->size){
        return 1;
    }
    return fmd_string_comp(c1->str, c2->str);
}
void *cs_copy(cs_t *c1) {
    cs_t *copy = calloc(1, sizeof(cs_t));
    copy->str = fmd_string_copy(c1->str);
    copy->vertices = fmd_vector_copy(c1->vertices);
    return copy;
}

static fmd_fstruct_t fstruct_cs = {
        (fcomp) cs_comp,
        (fhash) cs_hash,
        (ffree) cs_free,
        (fcopy) cs_copy
};

void flatten_constraints_helper(void *key, void* value, void *p) {
    fmd_vector_t *constraints = (fmd_vector_t*)p;
    cs_t *constraint = calloc(1, sizeof(cs_t));
    constraint->str = key;
    constraint->vertices = value;
    fmd_vector_append(constraints, constraint);
}

fmd_vector_t *flatten_constraints(fmd_table_t *constraint_sets) {
    fmd_vector_t *constraints;
    fmd_vector_init(&constraints, FMD_VECTOR_INIT_SIZE, &fstruct_cs);
    fmd_table_traverse(constraint_sets, constraints, flatten_constraints_helper);
    fmd_vector_sort(constraints);
    return constraints;
}

int main() {
    srand(1234567);
    int_t depth = 41;
    fmd_graph_t *graph = generate_random_graph(100, 150, 5, 6);
    fmd_table_t *c_table = fmd_graph_extract_constraint_sets(graph, depth);
    fmd_vector_t *constraints = flatten_constraints(c_table);

    //for(int_t i = 0; i < constraints->size; i++) {
    //    cs_t *constraint = constraints->data[i];
    //    printf("%s: [ ", constraint->str->seq);
    //    for(int j = 0; j < constraint->vertices->size; j++) {
    //        printf("%ld ", (int_t)constraint->vertices->data[j]);
    //    }
    //    printf("]\n");
    //}

    int_t cardinality[depth];
    int_t count[depth];
    for(int_t i = 0; i < depth; i++) {
        cardinality[i] = count[i] = 0;
    }
    for(int_t i = 0; i < constraints->size; i++) {
        cs_t *constraint = constraints->data[i];
        cardinality[constraint->str->size-1] += constraint->vertices->size;
        ++count[constraint->str->size-1];
    }
    for(int_t i = 0; i < depth; i++) {
        printf("%d\t%d\t%d\t%4.4lf\n", i+1, cardinality[i], count[i], cardinality[i] / (double)count[i]);
    }



    return 0;
}