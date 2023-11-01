#include "gin_fmi.h"
#include <time.h>



int main(int argc, char *argv[]) {

    char *fname = argv[1];
    char iname[256];
    memset(iname, 0 , 256);

    strcat(iname, fname);
    strcat(iname, ".fmi");

    FILE *idx = fopen(iname,"r");
    gin_fmi_t *fmi;
    if(!idx) {
        char *buffer;
        FILE *file = fopen(fname, "r");
        fprintf(stdout, "Index not found in path %s, constructing...\n", iname);
        // parse contents of fname
        // Seek to the end of the file to determine its size
        fseek(file, 0, SEEK_END);
        long file_size = ftell(file);
        fseek(file, 0, SEEK_SET);  // Reset file pointer to the beginning

        // Allocate memory to store the file contents
        buffer = (char *) malloc(file_size + 1);
        if (!buffer) {
            perror("Failed to allocate memory");
            fclose(file);
            return 1;
        }

        // Read the entire file into the buffer
        fread(buffer, 1, file_size, file);
        buffer[file_size] = '\0';  // Null-terminate if you want to treat it as a string
        gin_string_t *str;
        gin_string_init_cstr(&str,buffer);
        free(buffer);
        gin_fmi_init(&fmi,str,32,32);
        gin_string_free(str);
    }
    // time queries
    struct timespec t1,t2;
    double time = 0;
    int no_queries = 0;
    char line[16384];
    if (fgets(line, sizeof(line), stdin)) {
        // Remove newline character if it exists
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n') {
            line[len - 1] = '\0';
        }
        gin_string_t *query;
        gin_string_init_cstr(&query, line);
        clock_gettime(CLOCK_REALTIME, &t1);
        gin_fmi_query_count(fmi, query);
        clock_gettime(CLOCK_REALTIME, &t2);
        time += t2.tv_sec - t1.tv_sec + (t2.tv_nsec - t1.tv_nsec) * 1e-9;
        no_queries++;
    }
    fprintf(stdout, "Queries per second: %lf\n", no_queries / time);
    return 0;
}