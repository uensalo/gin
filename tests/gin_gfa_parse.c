#include "gin_gin.h"

#include <stdio.h>
#include <rgfa_parser.h>


int main(int argc, char *args[]) {
    printf("%s\n", args[1]);
    FILE *fp = fopen(args[1], "r");
    if(!fp) {
        fprintf(stderr,"Failed to open %s, quitting\n", args[1]);
        exit(-1);
    }
    printf("Parsing RGFA file under %s\n", args[1]);
    rgfa_t *rgfa = rgfa_parse(fp);
    printf("Number of segments: %d\n", rgfa->num_slines);
    printf("Number of edges:    %d\n", rgfa->num_llines);
    int_t no_chars = 0;
    for(int_t i = 0; i < rgfa->num_slines; i++) {
        no_chars += strlen(rgfa->slines[i].seq);
    }
    printf("total parsed chars: %ld\n", no_chars);

    gin_graph_t *g = rgfa_to_gin_graph(rgfa);
    printf("Number of segments in the graph: %d\n", g->vertex_list->size);
    printf("Number of edges in the graph:    %d\n", g->no_edges);
    no_chars = 0;
    for(int_t i = 0; i < g->vertex_list->size; i++) {
        gin_vertex_t* v;
        gin_table_lookup(g->vertices, (void*)i, &v);
        assert(strlen(rgfa->slines[i].seq) == v->label->size);
        no_chars += v->label->size;
    }
    printf("total chars: %ld\n", no_chars);

    //gin_gin_t *gin;
    //gin_gin_init(&gin,g,NULL,GIN_GIN_DEFAULT_c_0,GIN_GIN_DEFAULT_c_1,256,256);

    rgfa_free(rgfa);
    gin_graph_free(g);
    return 0;
}