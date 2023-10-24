/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * rgfa_parser.h is part of gin
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
#ifndef GIN_RGFA_PARSER_H
#define GIN_RGFA_PARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gin_graph.h>

typedef struct {
    char* segId;
    char* seq;
    char* sn;
    int so;
    int sr;
} rgfa_sline_t;

typedef struct {
    char* segId1;
    char strand1;
    char* segId2;
    char strand2;
    char* cigar;
} rgfa_lline_t;

typedef struct {
    rgfa_sline_t* slines;
    size_t num_slines;
    rgfa_lline_t* llines;
    size_t num_llines;
} rgfa_t;

static rgfa_t* rgfa_parse(FILE* file) {
    rgfa_t* rgfa = (rgfa_t*)calloc(1, sizeof(rgfa_t));
    long unsigned int slines_size = 8;
    long unsigned int llines_size = 8;
    rgfa->slines = calloc(slines_size, sizeof(rgfa_sline_t));
    rgfa->llines = calloc(llines_size, sizeof(rgfa_lline_t));

    size_t line_size = 4096;
    size_t buffer_size = line_size;
    char* line = calloc(line_size, sizeof(char));

    char *p, *end;

    while (1) {
        p = fgets(line, (int)line_size, file);

        // Stop if we have reached the end of the file.
        if (!p || feof(file)) {
            break;
        }

        // Expand the line buffer if needed.
        while (line[strlen(line) - 1] != '\n') {
            buffer_size *= 2;
            line = realloc(line, buffer_size);
            fgets(line + line_size - 1, (int)line_size + 1, file);
            line_size = buffer_size;
        }

        p = line + 2;  // Skip the first two characters ('S' or 'L' and '\t')

        if (line[0] == 'S') {
            rgfa_sline_t sline;
            sline.segId = NULL;
            sline.seq = NULL;
            sline.sn = NULL;

            // Parse segId
            end = strchr(p, '\t');
            sline.segId = strndup(p, end - p);
            p = end + 1;

            // Parse seq
            end = strchr(p, '\t');
            sline.seq = strndup(p, end - p);
            p = end + 1;

            // Parse tags
            while ((end = strchr(p, '\t')) != NULL) {
                if (p[0] == 'S' && p[1] == 'N') {
                    sline.sn = strndup(p + 5, end - p - 5);
                } else if (p[0] == 'S' && p[1] == 'O') {
                    sscanf(p, "SO:i:%d", &sline.so);
                } else if (p[0] == 'S' && p[1] == 'R') {
                    sscanf(p, "SR:i:%d", &sline.sr);
                }
                p = end + 1;
            }
            ++rgfa->num_slines;
            if(rgfa->num_slines >= slines_size) {
                rgfa->slines = (rgfa_sline_t *) realloc(rgfa->slines, (slines_size*=2) * sizeof(rgfa_sline_t));
            }
            rgfa->slines[rgfa->num_slines - 1] = sline;
        } else if (line[0] == 'L') {
            rgfa_lline_t lline;
            lline.segId1 = NULL;
            lline.segId2 = NULL;
            lline.cigar = NULL;

            // Parse segId1
            end = strchr(p, '\t');
            lline.segId1 = strndup(p, end - p);
            p = end + 1;

            // Parse strand1
            lline.strand1 = *p;
            p += 2;  // Skip the character and the next '\t'

            // Parse segId2
            end = strchr(p, '\t');
            lline.segId2 = strndup(p, end - p);
            p = end + 1;

            // Parse strand2
            lline.strand2 = *p;
            p += 2;  // Skip the character and the next '\t'

            // Parse cigar
            end = strchr(p, '\n');
            lline.cigar = strndup(p, end - p);

            ++rgfa->num_llines;
            if(rgfa->num_llines >= llines_size) {
                rgfa->llines = (rgfa_lline_t *) realloc(rgfa->llines, (llines_size*=2) * sizeof(rgfa_lline_t));
            }
            rgfa->llines[rgfa->num_llines - 1] = lline;
        }
    }
    free(line);
    rgfa->slines = (rgfa_sline_t *) realloc(rgfa->slines, rgfa->num_slines * sizeof(rgfa_sline_t));
    rgfa->llines = (rgfa_lline_t *) realloc(rgfa->llines, rgfa->num_llines * sizeof(rgfa_lline_t));
    return rgfa;
}

static void rgfa_free(rgfa_t * rgfa) {
    if (rgfa->slines) {
        for (size_t i = 0; i < rgfa->num_slines; ++i) {
            free(rgfa->slines[i].segId);
            free(rgfa->slines[i].seq);
            free(rgfa->slines[i].sn);
        }
        free(rgfa->slines);
    }
    if (rgfa->llines) {
        for (size_t i = 0; i < rgfa->num_llines; ++i) {
            free(rgfa->llines[i].segId1);
            free(rgfa->llines[i].segId2);
            free(rgfa->llines[i].cigar);
        }
        free(rgfa->llines);
    }
    free(rgfa);
}

static gin_graph_t *rgfa_to_gin_graph(rgfa_t *rgfa) {
    gin_graph_t* graph;
    gin_graph_init(&graph);

    for (int_t i = 0; i < rgfa->num_slines; i++) {
        rgfa_sline_t sline = rgfa->slines[i];
        gin_string_t* label;
        gin_string_init_cstr(&label, sline.seq);
        gin_graph_insert_vertex(graph, (vid_t)(atoi(sline.segId+1)-1), label);
    }

    for (int_t i = 0; i < rgfa->num_llines; i++) {
        rgfa_lline_t lline = rgfa->llines[i];
        int_t source_id = (vid_t)(atoi(lline.segId1+1)-1);
        int_t destination_id = (vid_t)(atoi(lline.segId2+1)-1);
        gin_graph_insert_edge(graph, source_id, destination_id);
    }

    return graph;
}

#endif //GIN_RGFA_PARSER_H
