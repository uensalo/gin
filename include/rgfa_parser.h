#ifndef FMD_RGFA_PARSER_H
#define FMD_RGFA_PARSER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fmd_graph.h>

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

rgfa_t* rgfa_parse(FILE* file) {
    rgfa_t* rgfa = (rgfa_t*)calloc(1, sizeof(rgfa_t));
    long unsigned int slines_size = 8;
    long unsigned int llines_size = 8;
    rgfa->slines = calloc(slines_size, sizeof(rgfa_sline_t));
    rgfa->llines = calloc(llines_size, sizeof(rgfa_lline_t));

    size_t line_size = 4096;
    size_t buffer_size = line_size;
    char* line = malloc(line_size);

    while (1) {
        line[0] = '\0';
        fgets(line, line_size, file);
        if (feof(file)) {
            break;
        }

        while (line[strlen(line) - 1] != '\n') {
            buffer_size *= 2;
            line = realloc(line, buffer_size);
            fgets(line + line_size - 1, line_size + 1, file);
            line_size = buffer_size;
        }

        if (line[0] == 'S') {
            rgfa_sline_t sline;
            sline.segId = NULL;
            sline.seq = NULL;
            sline.sn = NULL;

            char* tags = NULL;
            sscanf(line, "S\t%ms\t%ms\t%ms", &sline.segId, &sline.seq, &tags);

            char* tag = strtok(tags, "\t");
            while (tag != NULL) {
                if (tag[0] == 'S' && tag[1] == 'N') {
                    sscanf(tag, "SN:Z:%ms", &sline.sn);
                } else if (tag[0] == 'S' && tag[1] == 'O') {
                    sscanf(tag, "SO:i:%d", &sline.so);
                } else if (tag[0] == 'S' && tag[1] == 'R') {
                    sscanf(tag, "SR:i:%d", &sline.sr);
                }
                tag = strtok(NULL, "\t");
            }
            free(tags);
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

            sscanf(line, "L\t%ms\t%c\t%ms\t%c\t%ms", &lline.segId1, &lline.strand1, &lline.segId2, &lline.strand2, &lline.cigar);
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

void rgfa_free(rgfa_t * rgfa) {
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

#endif //FMD_RGFA_PARSER_H
