#include "gin_graph.h"
#include "gin_constraint_set.h"

gin_graph_t *g1() {
    gin_string_t *l0,*l1,*l2,*l3,*l4,*l5;
    gin_string_init_cstr(&l0, "A");
    gin_string_init_cstr(&l1, "C");
    gin_string_init_cstr(&l2, "G");
    gin_string_init_cstr(&l3, "T");
    gin_string_init_cstr(&l4, "A");
    gin_string_init_cstr(&l5, "C");

    gin_graph_t *g;
    gin_graph_init(&g);
    gin_graph_insert_vertex(g, 0, l0);
    gin_graph_insert_vertex(g, 1, l1);
    gin_graph_insert_vertex(g, 2, l2);
    gin_graph_insert_vertex(g, 3, l3);
    gin_graph_insert_vertex(g, 4, l4);
    gin_graph_insert_vertex(g, 5, l5);

    for(int_t i = 0; i < 6; i++) {
        for(int_t j = 0; j < i; j++) {
            gin_graph_insert_edge(g, i, j);
            gin_graph_insert_edge(g, j, i);
            printf("%d,%d\n",i,j);
        }
    }
    return g;
}


int main() {

    gin_graph_t *g = g1();
    gin_vector_t *constraints;
    gin_constraint_set_enumerate(&constraints, g, 4, true);

    for(int_t i = 0; i < constraints->size; i++) {
        gin_constraint_set_t *s = constraints->data[i];
        printf("%s : [ ", s->str->seq);
        for(int_t j = 0; j < s->vertices->size; j++) {
            printf("%d ", s->vertices->data[j]);
        }
        printf("]\n");
    }

}