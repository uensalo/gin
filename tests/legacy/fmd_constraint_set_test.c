#include "fmd_graph.h"
#include "fmd_constraint_set.h"

fmd_graph_t *g1() {
    fmd_string_t *l0,*l1,*l2,*l3,*l4,*l5;
    fmd_string_init_cstr(&l0, "A");
    fmd_string_init_cstr(&l1, "C");
    fmd_string_init_cstr(&l2, "G");
    fmd_string_init_cstr(&l3, "T");
    fmd_string_init_cstr(&l4, "A");
    fmd_string_init_cstr(&l5, "C");

    fmd_graph_t *g;
    fmd_graph_init(&g);
    fmd_graph_insert_vertex(g, 0, l0);
    fmd_graph_insert_vertex(g, 1, l1);
    fmd_graph_insert_vertex(g, 2, l2);
    fmd_graph_insert_vertex(g, 3, l3);
    fmd_graph_insert_vertex(g, 4, l4);
    fmd_graph_insert_vertex(g, 5, l5);

    for(int_t i = 0; i < 6; i++) {
        for(int_t j = 0; j < i; j++) {
            fmd_graph_insert_edge(g, i, j);
            fmd_graph_insert_edge(g, j, i);
            printf("%d,%d\n",i,j);
        }
    }
    return g;
}


int main() {

    fmd_graph_t *g = g1();
    fmd_vector_t *constraints;
    fmd_constraint_set_enumerate(&constraints, g, 4, true);

    for(int_t i = 0; i < constraints->size; i++) {
        fmd_constraint_set_t *s = constraints->data[i];
        printf("%s : [ ", s->str->seq);
        for(int_t j = 0; j < s->vertices->size; j++) {
            printf("%d ", s->vertices->data[j]);
        }
        printf("]\n");
    }

}