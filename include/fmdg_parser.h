#ifndef FMD_FMDG_PARSER_H
#define FMD_FMDG_PARSER_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fmd_graph.h"

static fmd_graph_t *fmdg_parse(FILE* file) {
    char *line = NULL;
    size_t len = 0;

    fmd_graph_t *graph;
    fmd_graph_init(&graph);

    while(getline(&line, &len, file) != -1) {
        char *start = line;
        while(*start == ' ' || *start == '\t') ++start;
        char *end = start + strlen(start) - 1;
        while(end > start && (*end == ' ' || *end == '\t' || *end == '\n')) --end;
        *(end + 1) = 0;

        // Skip empty lines
        if(*start == '\0') continue;

        char *type = strtok(line, "\t");
        if(type == NULL) continue;

        if(*type == 'V') {
            char_t *vid_str = strtok(NULL, "\t");
            char_t *vlabel_cstr = strtok(NULL, "\t");
            if(vid_str != NULL && vlabel_cstr != NULL) {
                int_t vid = strtoll(vid_str, NULL, 10);
                fmd_string_t *label;
                fmd_string_init_cstr(&label, vlabel_cstr);
                fmd_graph_insert_vertex(graph, vid, label);
            }
        }
        else if(*type == 'E') {
            char_t *vid1_str = strtok(NULL, "\t");
            char_t *vid2_str = strtok(NULL, "\t");
            if(vid1_str != NULL && vid2_str != NULL) {
                int_t vid1 = strtoll(vid1_str, NULL, 10);
                int_t vid2 = strtoll(vid2_str, NULL, 10);
                fmd_graph_insert_edge(graph, vid1, vid2);
            }
        }
    }
    return graph;
}



#endif //FMD_FMDG_PARSER_H
