#include "fmd_fmd.h"

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
    rgfa_free(rgfa);
    return 0;
}