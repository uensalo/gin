/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * ging_parser.h is part of gin
 *
 * gin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GIN_GING_PARSER_H
#define GIN_GING_PARSER_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gin_graph.h"

static gin_graph_t *ging_parse(FILE* file) {
    char *line = NULL;
    size_t len = 0;

    gin_graph_t *graph;
    gin_graph_init(&graph);

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
                gin_string_t *label;
                gin_string_init_cstr(&label, vlabel_cstr);
                gin_graph_insert_vertex(graph, vid, label);
            }
        }
        else if(*type == 'E') {
            char_t *vid1_str = strtok(NULL, "\t");
            char_t *vid2_str = strtok(NULL, "\t");
            if(vid1_str != NULL && vid2_str != NULL) {
                int_t vid1 = strtoll(vid1_str, NULL, 10);
                int_t vid2 = strtoll(vid2_str, NULL, 10);
                gin_graph_insert_edge(graph, vid1, vid2);
            }
        }
    }
    return graph;
}

static void *ging_write(FILE* file, gin_graph_t *graph) {
    for (int_t i = 0; i < graph->vertices->size; i++) {
        gin_vertex_t *vertex; void* _;
        gin_table_lookup(graph->vertices, (void*)i, &_);
        vertex = (gin_vertex_t*)_;
        fprintf(file, "V\t%lld\t%s\n", i, vertex->label->seq);
    }
    for (int_t i = 0; i < graph->vertices->size; i++) {
        gin_vector_t *outgoing; void* _;
        gin_table_lookup(graph->outgoing_neighbors, (void*)i, &_);
        outgoing = (gin_vector_t*)_;
        for (int_t j = 0; j < outgoing->size; j++) {
            vid_t vid = (vid_t)outgoing->data[j];
            fprintf(file, "E\t%lld\t%lld\n", i, vid);
        }
    }
}



#endif //GIN_GING_PARSER_H
