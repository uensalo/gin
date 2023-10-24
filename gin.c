/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin.c is part of gin
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
#include "gin_gin.h"

#include "gin_annealing.h"
#include <getopt.h>
#include <time.h>
#include "rgfa_parser.h"
#include "permutation_parser.h"
#include "ging_parser.h"
#include "gin_encoded_graph.h"
#ifdef GIN_OMP
#include <omp.h>
#endif

#define GIN_MAIN_ISA_SAMPLE_RATE_DEFAULT 32
#define GIN_MAIN_RANK_SAMPLE_RATE_DEFAULT 32
#define to_sec(t1,t2) ((double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9)
#define GIN_MAIN_BUF_READ_SIZE 1024

#define GIN_MAIN_QUERY_BUF_LEN 65536
#define GIN_MAIN_QUERY_EXIT_PROMPT "exit();"

char *gin_version = "1.0";
char* gin_mode_names[] = {"index","query","decode","permutation","utils","help"};

typedef enum gin_mode_ {
    gin_mode_index=0,
    gin_mode_query=1,
    gin_mode_decode=2,
    gin_mode_permutation=3,
    gin_mode_utils=4,
    gin_mode_help=5,
    gin_mode_no_modes=6
} gin_mode_t;

char* gin_query_mode_names[] = {"find","cache"};

typedef enum gin_query_mode_ {
    gin_query_mode_find=0,
    gin_query_mode_cache=1,
    gin_query_mode_no_modes=2,
} gin_query_mode_t;

char* gin_decode_mode_names[] = {"encode", "walks"};

typedef enum gin_decode_mode_ {
    gin_decode_mode_encode=0,
    gin_decode_mode_walks=1,
    gin_decode_mode_no_modes=2,
} gin_decode_mode_t;

char* gin_convert_mode_names[] = {"rgfa2ging", "fastq2query", "spectrum", "find"};

typedef enum gin_utils_mode_ {
    gin_utils_mode_rgfa2ging=0,
    gin_utils_mode_fastq2query=1,
    gin_utils_mode_spectrum=2,
    gin_utils_mode_find=3,
    gin_utils_mode_no_modes=4,
} gin_utils_mode_t;

int gin_main_index(int argc, char **argv);
int gin_main_query(int argc, char **argv, gin_query_mode_t mode);
int gin_main_decode(int argc, char **argv, gin_decode_mode_t mode);
int gin_main_permutation(int argc, char **argv);
int gin_main_utils(int argc, char **argv, gin_utils_mode_t mode);
int gin_main_help(gin_mode_t progmode, char *progname);

gin_mode_t gin_string_to_mode(char *str) {
    if(!str) return gin_mode_no_modes;
    gin_mode_t progmode = gin_mode_no_modes;
    for(int i = 0; i < gin_mode_no_modes; i++) {
        if(strcmp(str,gin_mode_names[i]) == 0) {
            progmode = (gin_mode_t)i;
            break;
        }
    }
    return progmode;
}

gin_query_mode_t gin_string_to_qmode(char *str) {
    if(!str) return gin_query_mode_no_modes;
    gin_query_mode_t progmode = gin_query_mode_no_modes;
    for(int i = 0; i < gin_query_mode_no_modes; i++) {
        if(strcmp(str,gin_query_mode_names[i]) == 0) {
            progmode = (gin_query_mode_t)i;
            break;
        }
    }
    return progmode;
}

gin_decode_mode_t gin_string_to_dmode(char *str) {
    if(!str) return gin_decode_mode_no_modes;
    gin_decode_mode_t progmode = gin_decode_mode_no_modes;
    for(int i = 0; i < gin_decode_mode_no_modes; i++) {
        if(strcmp(str,gin_decode_mode_names[i]) == 0) {
            progmode = (gin_decode_mode_t)i;
            break;
        }
    }
    return progmode;
}

gin_utils_mode_t gin_string_to_cmode(char *str) {
    if(!str) return gin_utils_mode_no_modes;
    gin_utils_mode_t progmode = gin_utils_mode_no_modes;
    for(int i = 0; i < gin_utils_mode_no_modes; i++) {
        if(strcmp(str,gin_convert_mode_names[i]) == 0) {
            progmode = (gin_utils_mode_t)i;
            break;
        }
    }
    return progmode;
}


int main(int argc, char *argv[]) {
    gin_mode_t progmode = gin_string_to_mode(argv[1]);
    char **argvp = argv+1;
    int argcp = argc-1;

    int return_code = -1;
    switch (progmode) {
        case gin_mode_index: {
            return_code = gin_main_index(argcp, argvp);
            break;
        }
        case gin_mode_query : {
            char *query_mode_name = argcp > 1 ? argvp[1] : NULL;
            gin_query_mode_t query_mode = gin_string_to_qmode(query_mode_name);
            return_code = gin_main_query(argcp, argvp, query_mode);
            break;
        }
        case gin_mode_decode : {
            char *query_mode_name = argcp > 1 ? argvp[1] : NULL;
            gin_decode_mode_t decode_mode = gin_string_to_dmode(query_mode_name);
            return_code = gin_main_decode(argcp, argvp, decode_mode);
            break;
        }
        case gin_mode_permutation : {
            return_code = gin_main_permutation(argcp, argvp);
            break;
        }
        case gin_mode_utils: {
            char *utils_mode_name = argcp > 1 ? argvp[1] : NULL;
            gin_utils_mode_t utils_mode = gin_string_to_cmode(utils_mode_name);
            return_code = gin_main_utils(argcp, argvp, utils_mode);
            break;
        }
        case gin_mode_help: {
            char *progname = argcp > 1 ? argvp[1] : NULL;
            gin_mode_t help_mode = gin_string_to_mode(progname);
            return_code = gin_main_help(help_mode, progname);
            break;
        }
        default: {
            fprintf(stderr, "[gin:error] Program %s not recognized. Please run gin help for more information. Quitting.\n", argv[1]);
            break;
        }
    }
    return return_code;
}

int gin_main_index(int argc, char **argv) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fperm_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fperm = NULL;
    int_t isa_rate = GIN_MAIN_ISA_SAMPLE_RATE_DEFAULT;
    int_t rank_rate = GIN_MAIN_RANK_SAMPLE_RATE_DEFAULT;
    bool parse_rgfa = false;
    bool verbose = false;

    struct timespec t1;
    struct timespec t2;
    double index_time = .0;
    double parse_time = .0;
    double write_time = .0;

    int return_code = 0;
    static struct option options[] = {
            {"input",                  required_argument, NULL, 'i'},
            {"rgfa",                   no_argument,       NULL, 'g'},
            {"output",                 required_argument, NULL, 'o'},
            {"permutation",            required_argument, NULL, 'p'},
            {"isa-sample-rate",        required_argument, NULL, 's'},
            {"rank-sample-rate",       required_argument, NULL, 'r'},
            {"verbose",                no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "i:go:p:s:r:v", options, &optindex)) != -1) {
        switch(c) {
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[gin:index] Input path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'g': {
                parse_rgfa = true;
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'p': {
                fperm_path = optarg;
                fperm = fopen(fperm_path, "r");
                if(!fperm) {
                    fprintf(stderr, "[gin:index] Permutation path %s could not be opened, quitting.\n", fperm_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 's': {
                isa_rate = strtoll(optarg, NULL, 10);
                break;
            }
            case 'r': {
                rank_rate = strtoll(optarg, NULL, 10);
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[gin:index] Option %s not recognized, please see gin help index for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    gin_graph_t *graph = NULL;
    gin_vector_t *permutation = NULL;
    gin_gin_t *gin = NULL;
    unsigned char *gin_buf;
    uint64_t gin_buf_size;

    /**************************************************************************
     * 1 - Parse the input graph the variable graph. By the end of this block,
     * the variable graph should be populated and the files must be closed.
     *************************************************************************/
    if(verbose) {
        fprintf(stderr, "[gin:index] Parsing input file %s\n", finput_path);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    if(parse_rgfa) {
        rgfa_t *rgfa = rgfa_parse(finput);
        if(finput != stdin) fclose(finput);
        if (!rgfa) {
            fprintf(stderr, "[gin:index] Failed to parse rGFA file under %s, quitting.\n", finput_path);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);
            if(fperm) fclose(fperm);
            return_code = -1;
            return return_code;
        }
        graph = rgfa_to_gin_graph(rgfa);
        rgfa_free(rgfa);

    } else {
        graph = ging_parse(finput);
        if(!graph) {
            if (finput_path) fclose(finput);
            if (foutput_path) fclose(foutput);
            if (fperm) fclose(fperm);
            fprintf(stderr, "[gin:index] Malformed ging file, quitting.\n");
            return_code = -1;
            return return_code;
        }
    }
    /**************************************************************************
     * 2 - Parse the permutation. By the end of this block, the permutation
     * file must be closed, and the variable permutation should be populated
     *************************************************************************/
    if(verbose) {
        if(fperm_path) fprintf(stderr, "[gin:index] Input parsed. Parsing permutation file %s\n", fperm_path);
        else fprintf(stderr, "[gin:index] Input parsed. Using the identity permutation.\n");
    }
    if(fperm) { // can not be provided through stdin; has to be from some file.
        permutation = permutation_parse(fperm);
        if(!permutation) {
            fprintf(stderr, "[gin:index] Failed to parse permutation file under %s, quitting.\n", fperm_path);
            return_code = -1;
            if(foutput != stdout) fclose(foutput);
            gin_graph_free(graph);
            return return_code;
        }
        if(permutation->size != graph->vertices->size) {
            fprintf(stderr, "[gin:index] Permutation cardinality does not match graph vertex size. Quitting.\n");
            return_code = -1;
            if(foutput_path) fclose(foutput);
            gin_vector_free(permutation);
            gin_graph_free(graph);
            return return_code;
        }
    }
    clock_gettime(CLOCK_REALTIME, &t2);
    parse_time = to_sec(t1,t2);
    /**************************************************************************
    * 3 - Construct the FM-Directory with the provided parameters
    *************************************************************************/
    if(verbose) { // print some statistics
        int_t no_vertices = graph->vertices->size;
        fprintf(stderr, "[gin:index] Constructing GIN on a graph with:\n");
        fprintf(stderr, "\tNumber of vertices: %lld\n", no_vertices);
        fprintf(stderr, "\tNumber of edges: %lld\n", graph->no_edges);
        int_t no_chars = 0;
        for(int_t i = 0; i < graph->vertex_list->size; i++) {
            void* _;
            gin_table_lookup(graph->vertices, (void*)i, &_);
            gin_vertex_t *v = (gin_vertex_t*)_;
            no_chars += v->label->size;
        }
        fprintf(stderr, "\tAverage in-degree: %.3lf\n", (double)(graph->no_edges) / (double)no_vertices);
        fprintf(stderr, "\tTotal length of vertex labels: %lld\n", no_chars);
        fprintf(stderr, "\tAverage label length per vertex: %.3lf\n", (double)no_chars / (double)no_vertices);
        int_t enc_sz = (no_chars + 2 * no_vertices + gin_ceil_log2(no_vertices) * no_vertices + 1);
        int_t sa_sz = (int_t)sizeof(int_t) * enc_sz;
        int_t graph_sz = no_chars + ((int_t)sizeof(gin_vertex_t) + 2*(int_t)sizeof(gin_vertex_t*)) * no_vertices + (int_t)sizeof(int_t) * 2 * graph->no_edges;
        int_t perm_sz = permutation ? (int_t)sizeof(vid_t) * permutation->size : 0;
        int_t emu = sa_sz + enc_sz + graph_sz + perm_sz;
        fprintf(stderr, "[gin:index] Estimated peak memory requirement: %lld bytes (%.3lf GB)\n", emu, (double)emu * 1e-9);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    gin_gin_init(&gin,
                 graph,
                 permutation,
                 GIN_GIN_DEFAULT_c_0,
                 GIN_GIN_DEFAULT_c_1,
                 rank_rate,
                 isa_rate);
    clock_gettime(CLOCK_REALTIME, &t2);
    index_time = to_sec(t1,t2);
    gin_graph_free(graph);
    gin_vector_free(permutation);
    /**************************************************************************
    * 4 - Convert the FM-Index to a buffer and write to disk
    *************************************************************************/
    clock_gettime(CLOCK_REALTIME, &t1);
    gin_gin_serialize_to_buffer(gin, &gin_buf, &gin_buf_size);
    gin_gin_free(gin);
    if(!gin_buf || !gin_buf_size) {
        fprintf(stderr, "[gin:index] Could not convert graph into a bitstream, please send a PR. Quitting.\n");
        return_code = -1;
        if(foutput != stdout) fclose(foutput);
        return return_code;
    }

    if(foutput_path) {
        foutput = fopen(foutput_path, "w");
        if(!foutput) {
            fprintf(stderr, "[gin:index] Output path %s could not be opened, quitting.\n", foutput_path);
            free(gin_buf);
            return_code = -1;
            return return_code;
        }
    }

    fwrite(gin_buf,sizeof(unsigned char),gin_buf_size,foutput);
    if(foutput_path) fclose(foutput);
    clock_gettime(CLOCK_REALTIME, &t2);
    write_time = to_sec(t1,t2);
    if(verbose) {
        fprintf(stderr, "[gin:index] Resulting index size: %lld bytes (%.3lf GB).\n", gin_buf_size, (double) gin_buf_size * 1e-9);
        fprintf(stderr, "[gin:index] Timings in seconds: \n");
        fprintf(stderr, "\t Parsing:  %lf\n", parse_time);
        fprintf(stderr, "\t Indexing: %lf\n", index_time);
        fprintf(stderr, "\t Writing:  %lf\n", write_time);
        fprintf(stderr, "[gin:index] Total relevant runtime: %lf seconds.\n", parse_time + index_time + write_time);
    }

    free(gin_buf);
    return return_code;
}

int gin_main_query(int argc, char **argv, gin_query_mode_t mode) {
    int return_code = 0;
    if(mode == gin_query_mode_no_modes) {
        fprintf(stderr, "[gin:query] Invalid query mode %s supplied. Please run gin help query for more information. Quitting...\n", argv[0]);
        return_code = -1;
        return return_code;
    }

    char *fref_path = NULL;
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fcache_path = NULL;

    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fref = NULL;
    FILE *fcache = NULL;

    int_t num_threads = 1;
    int_t batch_size = 8;
    int_t cache_depth = 0;
    int_t max_forks = -1;
    int_t max_matches = -1;
    bool decode = false;
    bool parse_fastq = false;
    bool verbose = false;

    struct timespec t1;
    struct timespec t2;

    // statistics for path enumeration
    int_t queries_processed = 0;
    int_t queries_decoded = 0;
    double index_parse_time = .0;
    double query_time = .0;
    double query_decoding_time = .0;
    double cache_build_time = .0;

    int_t no_matching_forks = 0;
    int_t no_missing_forks = 0;
    int_t no_matching_count = 0; // to compare with DFS
    int_t no_decoded_count = 0; // count the number of decodes

    int_t no_forks_advanced = 0;

    static struct option options[] = {
            {"reference",  required_argument, NULL, 'r'},
            {"input",      required_argument, NULL, 'i'},
            {"fastq",      no_argument,       NULL, 'f'},
            {"output",     required_argument, NULL, 'o'},
            {"cache-depth",required_argument, NULL, 'c'},
            {"cache"      ,required_argument, NULL, 'C'},
            {"max-forks",  required_argument, NULL, 'm'},
            {"max-matches",required_argument, NULL, 'M'},
            {"decode",     no_argument,       NULL, 'd'},
            {"batch-size", required_argument, NULL, 'b'},
            {"threads",    required_argument, NULL, 'j'},
            {"verbose",    no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "r:i:fo:c:C:m:M:db:j:v", options, &optindex)) != -1) {
        switch(c) {
            case 'r': {
                fref_path = optarg;
                break;
            }
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[gin:query] Input query path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'f': {
                parse_fastq = true;
                break;
            }
            case 'c': {
                cache_depth = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'C': {
                fcache_path = optarg;
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'd': {
                decode = true;
                break;
            }
            case 'm': {
                max_forks = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'M': {
                max_matches = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'b': {
                batch_size = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'j': {
                num_threads = (int_t)strtoull(optarg, NULL, 10);
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[gin:query] Option %s not recognized, please see gin help query for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    switch(mode) {
        case gin_query_mode_cache: {
            /**************************************************************************
            * 1 - Read the graph index from disk, gin_buf and gin_buf_size should be
            * populated and fref be closed by the end of this block
            **************************************************************************/
            fref = fopen(fref_path, "r");
            if(!fref) {
                fprintf(stderr, "[gin:query] Could not open index reference file %s, quitting.\n", fref_path);
                return_code = -1;
                return return_code;
            }

            uint64_t gin_buf_capacity = 65536;
            uint64_t gin_buf_size = 0;
            uint8_t *gin_buf = calloc(gin_buf_capacity, sizeof(uint8_t));

            uint64_t read_size;
            uint8_t *read_buf = calloc(GIN_MAIN_BUF_READ_SIZE, sizeof(uint8_t));
            while((read_size = fread(read_buf,sizeof(uint8_t), GIN_MAIN_BUF_READ_SIZE, fref))) {
                if(read_size + gin_buf_size > gin_buf_capacity) {
                    gin_buf = realloc(gin_buf, (gin_buf_capacity *= 2)*sizeof(uint8_t));
                    memset(gin_buf + (gin_buf_capacity/2), 0, sizeof(uint8_t)*(gin_buf_capacity/2));
                }
                memcpy(gin_buf+gin_buf_size,read_buf,read_size);
                gin_buf_size+=read_size;
            }
            gin_buf = realloc(gin_buf, gin_buf_size*sizeof(uint8_t));
            free(read_buf);
            fclose(fref);

            /**************************************************************************
            * 2 - Parse the bitstream into data structures, by the end of this block
            * gin should be populated and gin_buf should be freed
            **************************************************************************/
            if(verbose) {
                fprintf(stderr, "[gin:query] Parsing reference index into data structures.\n");
            }
            gin_gin_t *gin;
            clock_gettime(CLOCK_REALTIME, &t1);
            gin_gin_serialize_from_buffer(&gin, gin_buf, gin_buf_size);
            clock_gettime(CLOCK_REALTIME, &t2);
            index_parse_time += to_sec(t1,t2);

            if(!gin) {
                fprintf(stderr, "[gin:query] Error encountered while parsing reference index. Quitting.\n");
                return_code = -1;
                return return_code;
            }
            /**************************************************************************
            * 3 - Construct the cache on the index, convert it into a buffer and dump
            * the contents to the output file. By the end of this block, cache should
            * be written to the output and foutput should be closed if specified.
            **************************************************************************/
            if(verbose) {
                fprintf(stderr, "[gin:query] Constructing a cache of depth %lld.\n", cache_depth);
            }
            gin_gin_cache_t *cache; unsigned char *cache_buf; uint64_t cache_buf_size;
            clock_gettime(CLOCK_REALTIME, &t1);
            gin_gin_cache_init(&cache, gin, cache_depth);
            clock_gettime(CLOCK_REALTIME, &t2);
            cache_build_time += to_sec(t1,t2);

            if(verbose) {
                fprintf(stderr, "[gin:query] Cache contains %lld strings.\n", cache->no_entries);
                fprintf(stderr, "[gin:query] Dumping cache of size %.4lf MB to output.\n", gin_gin_cache_size(cache) / 1e6);
            }

            clock_gettime(CLOCK_REALTIME, &t1);
            gin_gin_cache_serialize_to_buffer(cache, &cache_buf, &cache_buf_size);
            clock_gettime(CLOCK_REALTIME, &t2);
            cache_build_time += to_sec(t1,t2);
            gin_gin_cache_free(cache);
            gin_gin_free(gin);
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:query] Can't open output file %s, quitting.\n", foutput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            fwrite(cache_buf,sizeof(unsigned char),cache_buf_size,foutput);
            free(cache_buf);
            break;
        }
        case gin_query_mode_find: {
            /**************************************************************************
            * 1 - Read the graph index from disk, gin_buf and gin_buf_size should be
            * populated and fref be closed by the end of this block
            **************************************************************************/
            fref = fopen(fref_path, "r");
            if(!fref) {
                fprintf(stderr, "[gin:query] Could not open index reference file %s, quitting.\n", fref_path);
                return_code = -1;
                return return_code;
            }

            uint64_t gin_buf_capacity = 65536;
            uint64_t gin_buf_size = 0;
            uint8_t *gin_buf = calloc(gin_buf_capacity, sizeof(uint8_t));

            uint64_t read_size;
            uint8_t *read_buf = calloc(GIN_MAIN_BUF_READ_SIZE, sizeof(uint8_t));
            while((read_size = fread(read_buf,sizeof(uint8_t), GIN_MAIN_BUF_READ_SIZE, fref))) {
                if(read_size + gin_buf_size > gin_buf_capacity) {
                    gin_buf = realloc(gin_buf, (gin_buf_capacity *= 2)*sizeof(uint8_t));
                    memset(gin_buf + (gin_buf_capacity/2), 0, sizeof(uint8_t)*(gin_buf_capacity/2));
                }
                memcpy(gin_buf+gin_buf_size,read_buf,read_size);
                gin_buf_size+=read_size;
            }
            gin_buf = realloc(gin_buf, gin_buf_size*sizeof(uint8_t));
            free(read_buf);
            fclose(fref);

            /**************************************************************************
            * 2 - Parse the bitstream into data structures, by the end of this block
            * gin, gin_dec, gin_cache should be populated and gin_buf should be freed
            **************************************************************************/
            if(verbose) {
                fprintf(stderr, "[gin:query] Parsing reference index into data structures.\n");
            }
            gin_gin_t *gin;
            gin_gin_cache_t *gin_cache = NULL;
            gin_gin_decoder_t *gin_dec;
            clock_gettime(CLOCK_REALTIME, &t1);
            gin_gin_serialize_from_buffer(&gin, gin_buf, gin_buf_size);
            clock_gettime(CLOCK_REALTIME, &t2);
            index_parse_time += to_sec(t1,t2);

            if(!gin) {
                fprintf(stderr, "[gin:query] Error encountered while parsing reference index. Quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(decode) {
                gin_gin_decoder_init(&gin_dec, gin);
            }
            if(fcache_path) {
                fcache = fopen(fcache_path, "r");
                if(!fcache) {
                    fprintf(stderr, "[gin:query] Error encountered while parsing cache. Quitting.\n");
                    return_code = -1;
                    return return_code;
                }
                if(verbose) {
                    fprintf(stderr, "[gin:query] Parsing cache from %s.\n", fcache_path);
                }
                uint64_t gin_cache_buf_capacity = 65536;
                uint64_t gin_cache_buf_size = 0;
                uint8_t *gin_cache_buf = calloc(gin_cache_buf_capacity, sizeof(uint8_t));

                uint64_t cache_read_size;
                uint8_t *cache_read_buf = calloc(GIN_MAIN_BUF_READ_SIZE, sizeof(uint8_t));
                clock_gettime(CLOCK_REALTIME, &t1);
                while((cache_read_size = fread(cache_read_buf,sizeof(uint8_t), GIN_MAIN_BUF_READ_SIZE, fcache))) {
                    if(cache_read_size + gin_cache_buf_size > gin_cache_buf_capacity) {
                        gin_cache_buf = realloc(gin_cache_buf, (gin_cache_buf_capacity *= 2)*sizeof(uint8_t));
                        memset(gin_cache_buf + (gin_cache_buf_capacity/2), 0, sizeof(uint8_t)*(gin_cache_buf_capacity/2));
                    }
                    memcpy(gin_cache_buf+gin_cache_buf_size,cache_read_buf,cache_read_size);
                    gin_cache_buf_size+=cache_read_size;
                }
                gin_cache_buf = realloc(gin_cache_buf, gin_cache_buf_size*sizeof(uint8_t));
                free(cache_read_buf);
                fclose(fcache);
                gin_gin_cache_serialize_from_buffer(&gin_cache, gin_cache_buf, gin_cache_buf_size);
                clock_gettime(CLOCK_REALTIME, &t2);
                cache_build_time += to_sec(t1,t2);
                cache_depth = gin_cache->depth;
                if(verbose) {
                    fprintf(stderr, "[gin:query] Cache of depth %lld parsed.\n", gin_cache->depth, (double)gin_gin_cache_size(gin_cache) / 1e6);
                }
            }
            /**************************************************************************
            * 3 - Parse queries from the input stream and query depending on the mode
            **************************************************************************/
            if(verbose) {
                fprintf(stderr, "[gin:query] Parsing and launching queries with %d threads and %d batch size.\n", (int)num_threads, (int)batch_size);
            }
            if(finput_path) {
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr,"[gin:query] Can't open input file %s, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:query] Can't open output file %s, quitting.\n", foutput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            if(parse_fastq) {
                // not yet implemented
                fprintf(stderr, "[gin:query] Fastq parsing mode not yet implemented. Quitting.\n");
                return_code = -1;
                gin_gin_free(gin);
                return return_code;
            } else {
                int_t i;
                char *buf = calloc(GIN_MAIN_QUERY_BUF_LEN, sizeof(char));
                typedef struct query_task_ {
                    gin_string_t *str;
                    uint64_t count;
                    gin_vector_t *exact_matches;
                    gin_vector_t *partial_matches;
                    gin_gin_stats_t *stats;
                } query_task_t;
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                query_task_t *tasks = calloc(batch_size, sizeof(query_task_t));
                bool exit_flag = false;

                while(true) {
                    i = 0;
                    while (i < batch_size && fgets(buf, GIN_MAIN_QUERY_BUF_LEN, finput)) {
                        int_t len = (int_t)strlen(buf);
                        if(buf[len-1] == '\n') buf[len-1] = 0; // get rid of the end line character
                        if(!strcmp(buf,GIN_MAIN_QUERY_EXIT_PROMPT)) {
                            exit_flag = true;
                            break;
                        }
                        gin_string_init_cstr(&tasks[i].str, buf);
                        queries_processed++;
                        i++;
                    }

                    clock_gettime(CLOCK_REALTIME, &t1);
                    #ifdef GIN_OMP
                    omp_set_num_threads((int)num_threads);
                    #endif
                    switch(mode) {
                        case gin_query_mode_find: {
                            #pragma omp parallel for default(none) shared(gin, gin_cache, i, tasks, num_threads, max_forks)
                            for(int_t k = 0; k < i; k++) {
                                gin_gin_query_find(gin, gin_cache, tasks[k].str, max_forks, &tasks[k].exact_matches, &tasks[k].partial_matches, &tasks[k].stats);
                            }
                            break;
                        }
                        default: {
                            break;
                        }
                    }
                    clock_gettime(CLOCK_REALTIME, &t2);
                    query_time += to_sec(t1,t2);

                    if(decode) {
                        switch(mode) {
                            case gin_query_mode_find: {
                                for(int_t j = 0; j < i; j++) {
                                    if(!tasks[j].str) continue;
                                    gin_vector_t *decoded_matches;
                                    no_matching_forks += tasks[j].exact_matches->size;
                                    no_missing_forks += tasks[j].partial_matches->size;
                                    clock_gettime(CLOCK_REALTIME, &t1);
                                    gin_gin_decoder_decode_ends(gin_dec,tasks[j].exact_matches, max_matches, &decoded_matches);
                                    clock_gettime(CLOCK_REALTIME, &t2);
                                    query_decoding_time += to_sec(t1, t2);
                                    if(verbose) {
                                        int_t matching_count = 0;
                                        for (int_t k = 0; k < tasks[j].exact_matches->size; k++) {
                                            gin_fork_node_t *fork = tasks[j].exact_matches->data[k];
                                            matching_count += fork->sa_hi - fork->sa_lo;
                                        }
                                        no_matching_count += matching_count;
                                        no_forks_advanced += tasks[j].stats->no_calls_to_advance_fork;
                                        fprintf(foutput, "%s:(c:%lld):(mf:%lld,pf:%lld):(a:%lld)\n",
                                                tasks[j].str->seq,
                                                matching_count,
                                                tasks[j].exact_matches->size,
                                                tasks[j].partial_matches->size,
                                                tasks[j].stats->no_calls_to_advance_fork);
                                    } else {
                                        for (int_t k = 0; k < tasks[j].exact_matches->size; k++) {
                                            gin_fork_node_t *fork = tasks[j].exact_matches->data[k];
                                            no_matching_count += fork->sa_hi - fork->sa_lo;
                                        }
                                        fprintf(foutput, "%s:\n", tasks[j].str->seq);
                                    }
                                    for(int_t k = 0; k < decoded_matches->size; k++) {
                                        gin_vector_t *decoded_matches_for_fork = decoded_matches->data[k];
                                        no_decoded_count += decoded_matches_for_fork->size;
                                        for(int_t l = 0; l < decoded_matches_for_fork->size; l++) {
                                            gin_decoded_match_t *decoded_match = decoded_matches_for_fork->data[l];
                                            fprintf(foutput, "\t(v:%lld,o:%lld)\n", decoded_match->vid, decoded_match->offset);
                                        }
                                    }
                                    if(!decoded_matches->size) {
                                        fprintf(foutput, "\t-\n");
                                    }
                                    gin_vector_free(decoded_matches);
                                    gin_vector_free(tasks[j].exact_matches);
                                    gin_vector_free(tasks[j].partial_matches);
                                    gin_string_free(tasks[j].str);
                                    free(tasks[j].stats);
                                }
                                break;
                            }
                            default: {
                                break;
                            }
                        }
                    }
                    else {
                        switch (mode) {
                            case gin_query_mode_find: {
                                for (int_t j = 0; j < i; j++) {
                                    if (!tasks[j].str) continue;
                                    int_t matching_count = 0;
                                    for (int_t k = 0; k < tasks[j].exact_matches->size; k++) {
                                        gin_fork_node_t *fork = tasks[j].exact_matches->data[k];
                                        matching_count += fork->sa_hi - fork->sa_lo;
                                    }
                                    no_matching_count += matching_count;
                                    no_forks_advanced += tasks[j].stats->no_calls_to_advance_fork;
                                    if(verbose) {
                                        fprintf(foutput, "%s:(c:%lld):(mf:%lld,pf:%lld):(a:%lld)\n",
                                                tasks[j].str->seq,
                                                matching_count,
                                                tasks[j].exact_matches->size,
                                                tasks[j].partial_matches->size,
                                                tasks[j].stats->no_calls_to_advance_fork);
                                    } else {
                                        fprintf(foutput, "%s:\n", tasks[j].str->seq);
                                    }
                                    no_matching_forks += tasks[j].exact_matches->size;
                                    no_missing_forks += tasks[j].partial_matches->size;
                                    for (int_t k = 0; k < tasks[j].exact_matches->size; k++) {
                                        gin_fork_node_t *fork = tasks[j].exact_matches->data[k];
                                        fprintf(foutput, "\t(%lld,%lld)\n", fork->sa_lo,fork->sa_hi);
                                    } if(!tasks[j].exact_matches->size) {
                                        fprintf(foutput, "\t-\n");
                                    }

                                    gin_vector_free(tasks[j].exact_matches);
                                    gin_vector_free(tasks[j].partial_matches);
                                    gin_string_free(tasks[j].str);
                                    free(tasks[j].stats);
                                    tasks[j].str = NULL;
                                }
                            }
                            default: {
                                break;
                            }
                        }
                    }
                    if(exit_flag)
                        break;
                }
                free(buf);
                free(tasks);
            }
            /**************************************************************************
            * 4 - Free all data structures and return
            **************************************************************************/
            fprintf(foutput,"%s\n", GIN_MAIN_QUERY_EXIT_PROMPT);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);

            if(verbose) {
                fprintf(stderr, "[gin:query] Params: Index file name (-r): %s\n", fref_path);
                if(finput_path) {
                    fprintf(stderr, "[gin:query] Params: Query file name (-i): %s\n", finput_path);
                }
                fprintf(stderr, "[gin:query] Params: Read batch size (-b): %lld\n", batch_size);
                fprintf(stderr, "[gin:query] Params: Threads (-j): %lld\n", num_threads);
                fprintf(stderr, "[gin:query] Params: Maximum forks tracked (-m): %lld\n", max_forks);
                fprintf(stderr, "[gin:query] Params: Maximum matches decoded (-M): %lld\n", max_matches);
                fprintf(stderr, "[gin:query] Params: Cache path: (-C): %s\n", fcache_path);
                fprintf(stderr, "[gin:query] Index: Index parse time (s): %lf\n", index_parse_time);
                if(fcache_path) {
                    fprintf(stderr, "[gin:query] Cache: Cache parse time in seconds: %lf\n",cache_build_time);
                    fprintf(stderr, "[gin:query] Cache: Cache size in memory (MB): %.4lf\n", (double)gin_gin_cache_size(gin_cache) / 1e6);
                }
                if(decode) {
                    fprintf(stderr, "[gin:query] Decode: Total matches decoded: %lld\n", no_decoded_count);
                    fprintf(stderr, "[gin:query] Decode: Total decoding time (s): %lf\n", query_decoding_time);
                    fprintf(stderr, "[gin:query] Decode: Matches decoded per second: %lf\n", (double)no_decoded_count / query_decoding_time);
                    fprintf(stderr, "[gin:query] Decode: Time per match decode (s): %lf\n", query_decoding_time / (double)no_decoded_count);
                }
                if(queries_processed) {
                    fprintf(stderr, "[gin:query] Find: Total queries processed: %lld\n",queries_processed);
                    fprintf(stderr, "[gin:query] Find: Total querying time (s): %lf\n",query_time);
                    fprintf(stderr, "[gin:query] Find: Queries per second: %lf\n",(double) queries_processed / (double) query_time);
                    fprintf(stderr, "[gin:query] Find: Time per query (s): %lf\n",(double) query_time / (double) queries_processed);
                }
                if((mode == gin_query_mode_find) && queries_processed && no_matching_forks) {
                    fprintf(stderr, "[gin:query] Find: Number of matching forks: %lld\n",no_matching_forks);
                    fprintf(stderr, "[gin:query] Find: Number of partial forks: %lld\n",no_missing_forks);
                    fprintf(stderr, "[gin:query] Find: Number of matches: %lld\n",no_matching_count);
                    fprintf(stderr, "[gin:query] Find: Number of fork advances: %lld\n",no_forks_advanced);
                    fprintf(stderr, "[gin:query] Find: Number of fork advances per query: %lf\n",(double) no_forks_advanced / (double) queries_processed);
                }
            }
            if(decode) {
                gin_gin_decoder_free(gin_dec);
            }
            if(fcache_path) {
                gin_gin_cache_free(gin_cache);
            }
            gin_gin_free(gin);
            break;
        }
        default: {
            fprintf(stderr, "[gin:query] Unrecognised mode for gin query, quitting. \n");
            return_code = -1;
            return return_code;
        }
    }

    return return_code;
}

int gin_main_decode(int argc, char **argv, gin_decode_mode_t mode) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fref_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fref = NULL;

    int_t batch_size = 8;
    int_t strings_processed = 0;
    int_t roots_processed = 0;
    int_t walks_processed = 0;

    double walk_process_time = 0;

    bool verbose = false;
    int_t num_threads = 1;
    static struct option options[] = {
            {"input",       required_argument, NULL, 'i'},
            {"output",      required_argument, NULL, 'o'},
            {"reference",   required_argument, NULL, 'r'},
            {"batch_size",  required_argument, NULL, 'b'},
            {"threads",     required_argument, NULL, 'j'},
            {"verbose",     no_argument,       NULL, 'v'},
    };
    int return_code = 0;
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "i:o:r:b:j:v", options, &optindex)) != -1) {
        switch(c) {
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[gin:decdode] Input path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'r': {
                fref_path = optarg;
                break;
            }
            case 'b': {
                batch_size = strtoll(optarg, NULL, 10);
                break;
            }
            case 'j': {
                num_threads = strtoll(optarg, NULL, 10);
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[gin:decode] Option %s not recognized, please see gin help permutation for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    gin_encoded_graph_t *encoded_graph = NULL;

    switch (mode) {
        case gin_decode_mode_encode: {
            gin_graph_t *graph = NULL;
            unsigned char *gin_encoded_graph_buf;
            uint64_t gin_encoded_graph_buf_size;
            /**************************************************************************
            * 1 - Parse the input graph the variable graph. By the end of this block,
            * the variable graph should be populated and the files must be closed.
            *************************************************************************/
            graph = ging_parse(finput);
            if(!graph) {
                if (finput_path) fclose(finput);
                if (foutput_path) fclose(foutput);
                fprintf(stderr, "[gin:decode] Malformed ging file, quitting.\n");
                return_code = -1;
                return return_code;
            }
            /**************************************************************************
            * 2 - log2 bit-encode the input graph. By the end of this block, the
            * variable graph should be freed and encoded graph should be populated.
            *************************************************************************/
            gin_encoded_graph_init(&encoded_graph, graph);
            if(!encoded_graph) {
                if (finput_path) fclose(finput);
                if (foutput_path) fclose(foutput);
                fprintf(stderr, "[gin:decode] Failed to encode ging file, quitting.\n");
                return_code = -1;
                return return_code;
            }
            gin_graph_free(graph);
            graph = NULL;
            /**************************************************************************
            * 3 - Convert the encoded graph to a payload and flush it on the disk
            *************************************************************************/
            gin_encoded_graph_serialize_to_buffer(encoded_graph, &gin_encoded_graph_buf, &gin_encoded_graph_buf_size);
            gin_encoded_graph_free(encoded_graph);
            encoded_graph = NULL;
            if(!gin_encoded_graph_buf || !gin_encoded_graph_buf_size) {
                fprintf(stderr, "[gin:decode] Could not convert graph into a bitstream, please send a PR. Quitting.\n");
                return_code = -1;
                if(foutput != stdout) fclose(foutput);
                return return_code;
            }

            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:index] Output path %s could not be opened, quitting.\n", foutput_path);
                    free(gin_encoded_graph_buf);
                    return_code = -1;
                    return return_code;
                }
            }

            fwrite(gin_encoded_graph_buf,sizeof(unsigned char),gin_encoded_graph_buf_size,foutput);
            if(foutput_path) fclose(foutput);
            free(gin_encoded_graph_buf);

            break;
        }
        case gin_decode_mode_walks : {
            /**************************************************************************
            * 1 - Read the encoded graph from disk, gin_buf and gin_buf_size should be
            * populated and fref be closed by the end of this block
            **************************************************************************/
            fref = fopen(fref_path, "r");
            if(!fref) {
                fprintf(stderr, "[gin:decode] Could not open encoded graph reference file %s, quitting.\n", fref_path);
                return_code = -1;
                return return_code;
            }

            uint64_t gin_buf_capacity = 65536;
            uint64_t gin_buf_size = 0;
            uint8_t *gin_buf = calloc(gin_buf_capacity, sizeof(uint8_t));

            uint64_t read_size;
            uint8_t *read_buf = calloc(GIN_MAIN_BUF_READ_SIZE, sizeof(uint8_t));
            while((read_size = fread(read_buf,sizeof(uint8_t), GIN_MAIN_BUF_READ_SIZE, fref))) {
                if(read_size + gin_buf_size > gin_buf_capacity) {
                    gin_buf = realloc(gin_buf, (gin_buf_capacity *= 2)*sizeof(uint8_t));
                    memset(gin_buf + (gin_buf_capacity/2), 0, sizeof(uint8_t)*(gin_buf_capacity/2));
                }
                memcpy(gin_buf+gin_buf_size,read_buf,read_size);
                gin_buf_size+=read_size;
            }
            gin_buf = realloc(gin_buf, gin_buf_size*sizeof(uint8_t));
            free(read_buf);
            fclose(fref);
            /**************************************************************************
            * 2 - Read the graph bitstream into data structures
            **************************************************************************/
            struct timespec t1,t2;
            clock_gettime(CLOCK_REALTIME, &t1);
            gin_encoded_graph_serialize_from_buffer(&encoded_graph, gin_buf, gin_buf_size);
            clock_gettime(CLOCK_REALTIME, &t2);
            double data_struct_read = to_sec(t1,t2);
            if(verbose) {
                fprintf(stderr, "[gin:decode] Encoded graph read into data structures in %lf seconds\n", data_struct_read);
            }
            /**************************************************************************
            * 3 - Parse queries from input and launch them as tasks
            **************************************************************************/
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:decode] Can't open output file %s, quitting.\n", foutput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            int_t i, j = 0;
            char *buf = calloc(GIN_MAIN_QUERY_BUF_LEN, sizeof(char));
            char *outbuf = calloc(GIN_MAIN_QUERY_BUF_LEN, sizeof(char));
            int_t outbuf_size = GIN_MAIN_QUERY_BUF_LEN;
            typedef struct walk_task_ {
                vid_t v;
                int_t o;
                gin_vector_t *walks;
            } walk_task_t;
            typedef struct walk_block_ {
                gin_string_t *str;
                gin_bs_t *encoded_str;
                walk_task_t *tasks;
                int_t no_tasks;
            } walk_block_t;
            #ifdef GIN_OMP
            omp_set_num_threads((int)num_threads);
            #endif
            walk_block_t *blocks = calloc(batch_size, sizeof(walk_block_t));
            for(int_t b = 0; b < batch_size; b++) {
                blocks[b].tasks = calloc(batch_size, sizeof(walk_task_t));
            }
            bool exit_flag = false;
            bool iter_first = true;
            bool carryover = false;
            gin_string_t *last_printed_string = NULL;

            while(true) {
                i = 0;
                j = 0;
                iter_first = true;
                while (i < batch_size && fgets(buf, GIN_MAIN_QUERY_BUF_LEN, finput)) {
                    int_t len = (int_t)strlen(buf);
                    if(buf[len-1] == '\n') buf[len-1] = 0; // get rid of the end line character
                    if(!strcmp(buf,GIN_MAIN_QUERY_EXIT_PROMPT)) {
                        exit_flag = true;
                        break;
                    }
                    if (buf[0] != '\t') { // end of previous block, start of a new block
                        // finalize previous block if necessary
                        if (!iter_first) {
                            // encode the string if necessary
                            if (blocks[j].no_tasks) {
                                gin_bs_t *encoded_query;
                                int_t bits_per_char = gin_ceil_log2(encoded_graph->alphabet_size);
                                int_t no_words_in_encoding = 1 + (blocks[j].str->size * bits_per_char - 1) / WORD_NUM_BITS;
                                gin_bs_init_reserve(&encoded_query, no_words_in_encoding);
                                word_t idx = 0;
                                for(int_t q = 0; q < blocks[j].str->size; q++) {
                                    gin_bs_write_word(encoded_query, idx, (word_t)encoded_graph->encoding_table[blocks[j].str->seq[q]], bits_per_char);
                                    idx += bits_per_char;
                                }
                                blocks[j].encoded_str = encoded_query;
                            }
                            ++j;
                        } else if (carryover) {
                            gin_string_free(blocks[j].str);
                            gin_bs_free(blocks[j].encoded_str);
                            blocks[j].encoded_str = NULL;
                        }
                        // parse the string for the new block
                        buf[len-2] = 0; // get rid of the colon
                        gin_string_init_cstr(&blocks[j].str, buf);
                        ++strings_processed;
                    } else { // v,o pair
                        vid_t v = -1;
                        int_t o = -1;
                        if (buf[0] == '\t'&& buf[1] == '-') {
                            ++roots_processed;
                            ++i;
                        }
                        else if (sscanf(buf, "\t(v:%llu,o:%lld)", &v, &o) == 2) {
                            blocks[j].tasks[blocks[j].no_tasks].v = v;
                            blocks[j].tasks[blocks[j].no_tasks].o = o;
                            blocks[j].tasks[blocks[j].no_tasks].walks = NULL;
                            ++blocks[j].no_tasks;
                            ++roots_processed;
                            ++i;
                        }
                    }
                    iter_first = false;
                }
                if (!iter_first) { // encode the last block
                    if(blocks[j].no_tasks) {
                        gin_bs_t *encoded_query;
                        int_t bits_per_char = gin_ceil_log2(encoded_graph->alphabet_size);
                        int_t no_words_in_encoding = 1 + (blocks[j].str->size * bits_per_char - 1) / WORD_NUM_BITS;
                        gin_bs_init_reserve(&encoded_query, no_words_in_encoding);
                        word_t idx = 0;
                        for (int_t q = 0; q < blocks[j].str->size; q++) {
                            gin_bs_write_word(encoded_query, idx,
                                              (word_t) encoded_graph->encoding_table[blocks[j].str->seq[q]],
                                              bits_per_char);
                            idx += bits_per_char;
                        }
                        blocks[j].encoded_str = encoded_query;
                    }
                    ++j;
                }


                clock_gettime(CLOCK_REALTIME, &t1);
                // flatten blocks to tasks
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                #pragma omp parallel default(none) shared(blocks, encoded_graph) firstprivate(j)
                {
                    #pragma omp single
                    {
                        for (int_t b = 0; b < j; b++) {
                            for (int_t t = 0; t < blocks[b].no_tasks; t++) {
                                #pragma omp task default(none) shared(blocks, encoded_graph) firstprivate(b, t)
                                {
                                    gin_encoded_graph_walk_string(encoded_graph,
                                                                  blocks[b].str,
                                                                  blocks[b].tasks[t].v,
                                                                  blocks[b].tasks[t].o,
                                                                  blocks[b].encoded_str,
                                                                  gin_encoded_graph_walk_extend_default,
                                                                  &blocks[b].tasks[t].walks);
                                }
                            }
                        }
                    }
                    #pragma omp taskwait
                }

                clock_gettime(CLOCK_REALTIME, &t2);
                walk_process_time += to_sec(t1,t2);
                for(int_t b = 0; b < j; b++) {
                    if(last_printed_string != blocks[b].str) {
                        fprintf(foutput, "%s:\n", blocks[b].str->seq);
                        last_printed_string = blocks[b].str;
                    }
                    if(blocks[b].no_tasks) {
                        for (int_t t = 0; t < blocks[b].no_tasks; t++) {
                            walks_processed += blocks[b].tasks[t].walks->size;
                            for (int_t l = 0; l < blocks[b].tasks[t].walks->size; l++) {
                                gin_walk_t *walk = blocks[b].tasks[t].walks->data[l];
                                int_t start_offset = walk->head->graph_lo;
                                int_t end_offset = walk->tail->graph_hi;
                                gin_walk_node_t *n = walk->head;
                                char tmp[256];
                                int_t outbuf_len = 0;
                                while (n != walk->dummy) {
                                    int_t s = sprintf(tmp, "%lld", n->vid);
                                    while (s + outbuf_len + 1 >= outbuf_size) {
                                        outbuf = realloc(outbuf, (outbuf_size *= 2));
                                    }
                                    memcpy(outbuf + outbuf_len, tmp, s);
                                    outbuf_len += s + 1;
                                    outbuf[outbuf_len - 1] = ':';
                                    n = n->next;
                                }
                                outbuf[outbuf_len - 1] = '\0';
                                fprintf(foutput, "\t(%lld,%lld);%s\n", start_offset, end_offset, outbuf);
                            }
                            gin_vector_free(blocks[b].tasks[t].walks);
                            blocks[b].tasks[t].walks = NULL;
                        }
                        blocks[b].no_tasks = 0;
                    } else {
                        fprintf(foutput, "\t-\n");
                    }
                    if(b != j-1) {
                        gin_string_free(blocks[b].str);
                        gin_bs_free(blocks[b].encoded_str);
                        blocks[b].str = NULL;
                        blocks[b].encoded_str = NULL;
                    }
                }

                // move the contents of the last block to the first index
                blocks[0].str = blocks[j-1].str;
                blocks[0].encoded_str = blocks[j-1].encoded_str;
                if(j > 1) { // self swap
                    blocks[j - 1].str = NULL;
                    blocks[j - 1].encoded_str = NULL;
                }
                carryover = true;
                if(exit_flag)
                    break;
            }

            free(buf);
            free(outbuf);
            for(int_t b = 0; b < batch_size; b++) {
                free(blocks[b].tasks);
            }
            free(blocks);
            gin_encoded_graph_free(encoded_graph);

            if(verbose) {
                fprintf(stderr, "[gin:decode] Walks: Number of threads: %lld\n",num_threads);
                fprintf(stderr, "[gin:decode] Walks: Number of strings processed: %lld\n",strings_processed);
                fprintf(stderr, "[gin:decode] Walks: Number of roots processed: %lld\n",roots_processed);
                fprintf(stderr, "[gin:decode] Walks: Number of matching walks: %lld\n",walks_processed);
                fprintf(stderr, "[gin:decode] Walks: Time elapsed assembling walks: %lf\n",walk_process_time);
                fprintf(stderr, "[gin:decode] Walks: Time per root: %lf\n",(double) walk_process_time / (double) roots_processed);
                fprintf(stderr, "[gin:decode] Walks: Roots per second: %lf\n",(double) roots_processed / (double) walk_process_time);
            }
            break;
        }
        default:
            break;
    }

    return return_code;
}

int gin_main_permutation(int argc, char **argv) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fperm_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fperm = NULL;

    bool parse_rgfa = false;
    bool multiple_vertex_span = true;

    int_t depth = 4;
    int_t time = 15;
    int_t update = 3;
    bool verbose = false;
    int_t num_threads = 1;

    struct timespec t1;
    struct timespec t2;

    double temperature = 1e6;
    double cooling_factor = 0.95;

    int return_code = 0;
    static struct option options[] = {
            {"input",       required_argument, NULL, 'i'},
            {"rgfa",        no_argument,       NULL, 'g'},
            {"output",      required_argument, NULL, 'o'},
            {"permutation", required_argument, NULL, 'p'},
            {"span-paths",  no_argument,       NULL, 's'},
            {"depth",       required_argument, NULL, 'd'},
            {"temperature", required_argument, NULL, 'e'},
            {"cooling",     required_argument, NULL, 'c'},
            {"time",        required_argument, NULL, 't'},
            {"update",      required_argument, NULL, 'u'},
            {"threads",     required_argument, NULL, 'j'},
            {"verbose",     no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "i:go:p:sd:e:c:t:u:j:v", options, &optindex)) != -1) {
        switch(c) {
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[gin:permutation] Input path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'g': {
                parse_rgfa = true;
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'p': {
                fperm_path = optarg;
                fperm = fopen(fperm_path, "r");
                if(!fperm) {
                    fprintf(stderr, "[gin:permutation] Permutation path %s could not be opened, quitting.\n", fperm_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 's' : {
                multiple_vertex_span = false;
                break;
            }
            case 'd': {
                depth = strtoll(optarg, NULL, 10);
                break;
            }
            case 'e': {
                temperature = strtod(optarg, NULL);
                break;
            }
            case 'c': {
                cooling_factor = strtod(optarg, NULL);
                break;
            }
            case 't': {
                time = strtoll(optarg, NULL, 10);
                break;
            }
            case 'u': {
                update = strtoll(optarg, NULL, 10);
                break;
            }
            case 'j': {
                num_threads = strtoll(optarg, NULL, 10);
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[gin:permutation] Option %s not recognized, please see gin help permutation for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    gin_graph_t *graph = NULL;
    gin_vector_t *permutation = NULL;
    gin_vector_t *constraints = NULL;
    gin_annealing_t *ann = NULL;
    gin_vector_t *optimized_permutation = NULL;

    /**************************************************************************
     * 1 - Parse the input graph the variable graph. By the end of this block,
     * the variable graph should be populated and the files must be closed.
     *************************************************************************/
    if(verbose) {
        fprintf(stderr, "[gin:permutation] Parsing input file %s\n", finput_path);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    if(parse_rgfa) {
        rgfa_t *rgfa = rgfa_parse(finput);
        if(finput != stdin) fclose(finput);
        if (!rgfa) {
            fprintf(stderr, "[gin:permutation] Failed to parse rGFA file under %s, quitting.\n", finput_path);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);
            if(fperm) fclose(fperm);
            return_code = -1;
            return return_code;
        }
        graph = rgfa_to_gin_graph(rgfa);
        rgfa_free(rgfa);

    } else {
        graph = ging_parse(finput);
        if(!graph) {
            if (finput_path) fclose(finput);
            if (foutput_path) fclose(foutput);
            if (fperm) fclose(fperm);
            fprintf(stderr, "[gin:permutation] Malformed ging file, quitting.\n");
            return_code = -1;
            return return_code;
        }
    }

    /**************************************************************************
     * 2 - Parse the initial input permutation if any
    **************************************************************************/
    if(fperm) { // can not be provided through stdin; has to be from some file.
        permutation = permutation_parse(fperm);
        if(!permutation) {
            fprintf(stderr, "[gin:permutation] Failed to parse permutation file under %s, quitting.\n", fperm_path);
            return_code = -1;
            if(foutput != stdout) fclose(foutput);
            gin_graph_free(graph);
            return return_code;
        }
        if(permutation->size != graph->vertices->size) {
            fprintf(stderr, "[gin:permutation] Permutation cardinality does not match graph vertex size. Quitting.\n");
            return_code = -1;
            if(foutput_path) fclose(foutput);
            gin_vector_free(permutation);
            gin_graph_free(graph);
            return return_code;
        }
    }

    /**************************************************************************
     * 3 - Enumerate constraint sets, by the end of this block constraints
     * must be populated with constraint sets.
     **************************************************************************/
    if(verbose) {
        fprintf(stderr, "[gin:permutation] Extracting constraint sets for each occurring prefix.\n");
        char *spanstr = multiple_vertex_span ? "True" : "False";
        fprintf(stderr, "[gin:permutation] Constraint sets can span multiple vertices: %s\n", spanstr);
    }
    gin_constraint_set_enumerate(&constraints, graph, depth, multiple_vertex_span);
    if(!constraints) {
        fprintf(stderr, "[gin:permutation] Something went wrong during constraint set extraction. Quitting.\n");
        gin_graph_free(graph);
        return_code = -1;
        return return_code;
    }

    /**************************************************************************
     * 4 - Configure the annealing optimizer and begin optimizing
     *************************************************************************/
#ifdef GIN_OMP
    omp_set_num_threads((int)num_threads);
#endif
    gin_annealing_configure(&ann, graph, constraints, permutation, temperature, 1, cooling_factor, 0);
    gin_graph_free(graph);
    gin_vector_free(constraints);
    if(verbose) {
        fprintf(stderr, "[gin:permutation] Optimization begins with initial cost %.3lf for depth=%lld for %lld seconds.\n", ann->cur_cost, depth, time);
    }
    int_t no_iterations = 1 + (time - 1) / update;
    for(int_t i = 0; i < no_iterations; i++) {
        gin_annealing_iterate_seconds(ann, (int_t)update);
        if(verbose) {
            fprintf(stderr, "\tIteration %d, best_cost = %.3lf, cur_cost = %.3lf\n", ann->cur_iter + 1, ann->best_cost_so_far, ann->cur_cost);
        }
    }

    /**************************************************************************
     * 5 - Cleanup and write to the output
     *************************************************************************/
    gin_annealing_get_permutation(ann, &optimized_permutation);
    //gin_annealing_free(ann);
    if(foutput_path) {
        foutput = fopen(foutput_path, "w");
        if(!foutput) {
            fprintf(stderr, "[gin:index] Output path %s could not be opened, quitting.\n", foutput_path);
            return_code = -1;
            return return_code;
        }
    }

    for(int_t i = 0; i < optimized_permutation->size; i++) {
        fprintf(foutput, "%lld\n", (int_t)optimized_permutation->data[i]);
    }

    if(foutput_path) fclose(foutput);
    gin_vector_free(permutation);
    gin_vector_free(optimized_permutation);

    return return_code;
}

int gin_main_utils(int argc, char **argv, gin_utils_mode_t mode) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fref_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fref = NULL;

    bool verbose = false;

    struct timespec t1;
    struct timespec t2;

    double parse_time = 0.0;
    double convert_time = 0.0;
    double write_time = 0.0;
    double utils_find_time = 0.0;

    int_t utils_queries_processed = 0;
    int_t utils_no_matches = 0;

    int_t kmer = 64;
    int_t num_threads = 1;
    int_t batch_size = 8;

    int return_code = 0;
    static struct option options[] = {
            {"input",   required_argument, NULL, 'i'},
            {"output",  required_argument, NULL, 'o'},
            {"kmer",    required_argument, NULL, 'k'},
            {"verbose", no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "i:o:r:k:j:v", options, &optindex)) != -1) {
        switch(c) {
            case 'i': {
                finput_path = optarg;
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'r': {
                fref_path = optarg;
                break;
            }
            case 'k': {
                kmer = strtoll(optarg, NULL, 10);
                break;
            }
            case 'j': {
                num_threads = strtoll(optarg, NULL, 10);
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[gin:utils] Option %s not recognized, please see gin help convert for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    gin_graph_t *graph = NULL;
    switch (mode) {
        case gin_utils_mode_rgfa2ging: {
            finput = fopen(finput_path, "r");
            if(!finput) {
                fprintf(stderr, "[gin:utils] Input path %s could not be opened, quitting.\n", finput_path);
                return_code = -1;
                return return_code;
            }
            // first parse the graph
            clock_gettime(CLOCK_REALTIME, &t1);
            if(verbose) {
                if(finput_path) fprintf(stderr, "[gin:utils] Parsing input file %s\n", finput_path);
                else fprintf(stderr, "[gin:utils] Parsing input from stdin.\n");
            }
            rgfa_t *rgfa = rgfa_parse(finput);
            clock_gettime(CLOCK_REALTIME, &t2);
            parse_time = to_sec(t1,t2);
            if(finput_path) fclose(finput);
            if(!rgfa) {
                fprintf(stderr, "[gin:utils] Failed to parse rgfa file quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(verbose) {
                if(finput_path) fprintf(stderr, "[gin:utils] Converting to ging format.\n");
            }
            clock_gettime(CLOCK_REALTIME, &t1);
            graph = rgfa_to_gin_graph(rgfa);
            clock_gettime(CLOCK_REALTIME, &t2);
            convert_time = to_sec(t1,t2);
            rgfa_free(rgfa);
            if(!graph) {
                fprintf(stderr, "[gin:utils] Failed to convert rgfa_t to gin_graph_t. Quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:utils] Failed to open output path %s, quitting\n", foutput_path);
                    gin_graph_free(graph);
                    return_code = -1;
                    return return_code;
                }
            }
            if(verbose) {
                if(foutput_path) fprintf(stderr, "[gin:utils] Writing converted graph to %s.\n", foutput_path);
                else fprintf(stderr, "[gin:utils] Writing converted graph to stdout.\n");
            }
            clock_gettime(CLOCK_REALTIME, &t1);
            ging_write(foutput, graph);
            clock_gettime(CLOCK_REALTIME, &t2);
            write_time = to_sec(t1,t2);
            if(foutput_path) fclose(foutput);
            if(verbose) {
                fprintf(stderr, "[gin:utils] Timings in seconds: \n");
                fprintf(stderr, "\t Parsing:    %.3lf\n", parse_time);
                fprintf(stderr, "\t Converting: %.3lf\n", convert_time);
                fprintf(stderr, "\t Writing:    %.3lf\n", write_time);
                fprintf(stderr, "[gin:utils] Total relevant runtime: %lf seconds.\n", parse_time + convert_time + write_time);
            }
            break;
        }
        case gin_utils_mode_fastq2query: {
            fprintf(stderr, "[gin:utils] fasta2query not implemented yet. Quitting.\n");
            return_code = -1;
            break;
        }
        case gin_utils_mode_spectrum: {
            finput = fopen(finput_path, "r");
            if(!finput) {
                fprintf(stderr, "[gin:utils] Input path %s could not be opened, quitting.\n", finput_path);
                return_code = -1;
                return return_code;
            }
            if(verbose) {
                fprintf(stderr, "[gin:utils] Parsing ging file %s\n", finput_path);
            }
            graph = ging_parse(finput);
            if(!graph) {
                if (finput_path) fclose(finput);
                fprintf(stderr, "[gin:utils] Malformed ging file, quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(finput_path) {
                fclose(finput);
            }
            if(verbose) {
                fprintf(stderr, "[gin:utils] ging file parsed, computing %lld-mer spectrum.\n", kmer);
            }
            gin_vector_t *spectrum;
            gin_table_t *spectrum_table;
            gin_graph_kmer_locations(graph, kmer, &spectrum, &spectrum_table);
            gin_graph_free(graph);
            if(verbose) {
                fprintf(stderr, "[gin:utils] Spectrum contains %lld %lld-mers. Sorting ocurrences.\n", spectrum_table->size, kmer);
            }
            gin_vector_sort(spectrum);
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:utils] Can't open output path %s\n", foutput_path);
                    gin_vector_free_disown(spectrum);
                    gin_table_free(spectrum_table);
                    return_code = -1;
                    return return_code;
                }
            }
            for(int_t i = 0; i < spectrum->size; i++) {
                gin_kmer_kv_t *kv = spectrum->data[i];
                fprintf(foutput, "%s:(v:%lld,o:%lld)\n", kv->str->seq, kv->vid, kv->offset);
            }
            if(foutput_path) {
                fclose(foutput);
            }
            gin_vector_free_disown(spectrum);
            gin_table_free(spectrum_table);
            break;
        }
        case gin_utils_mode_find: {
            if(verbose) {
                fprintf(stderr, "[gin:utils] Parsing ging file %s\n", fref_path);
            }
            fref = fopen(fref_path, "r");
            if(!fref) {
                fprintf(stderr, "[gin:utils] Reference path %s could not be opened, quitting.\n", fref);
                return_code = -1;
                return return_code;
            }
            clock_gettime(CLOCK_REALTIME, &t1);
            graph = ging_parse(fref);
            clock_gettime(CLOCK_REALTIME, &t2);
            parse_time += to_sec(t1,t2);
            if(!graph) {
                if (fref_path) fclose(fref);
                fprintf(stderr, "[gin:utils] Malformed ging file, quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(verbose) {
                fprintf(stderr, "[gin:utils] Graph file %s parsed.\n", fref_path);
            }
            if(verbose) {
                fprintf(stderr, "[gin:utils] Parsing and launching queries with %d threads.\n", (int)num_threads);
            }
            if(finput_path) {
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr,"[gin:utils] Can't open input file %s, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[gin:utils] Can't open output file %s, quitting.\n", foutput_path);
                    return_code = -1;
                    return return_code;
                }
            }
            int_t i;
            char *buf = calloc(GIN_MAIN_QUERY_BUF_LEN, sizeof(char*));
            typedef struct utils_find_task_ {
                gin_string_t *str;
                uint64_t count;
                gin_vector_t *matches;
            } utils_find_task_t;
            #ifdef GIN_OMP
            omp_set_num_threads((int)num_threads);
            #endif
            utils_find_task_t *tasks = calloc(batch_size, sizeof(utils_find_task_t));
            bool exit_flag = false;

            while(true) {
                i = 0;
                while (i < batch_size && fgets(buf, GIN_MAIN_QUERY_BUF_LEN, finput)) {
                    int_t len = (int_t)strlen(buf);
                    if(buf[len-1] == '\n') buf[len-1] = 0; // get rid of the end line character
                    if(!strcmp(buf,GIN_MAIN_QUERY_EXIT_PROMPT)) {
                        exit_flag = true;
                        break;
                    }
                    gin_string_init_cstr(&tasks[i].str, buf);
                    utils_queries_processed++;
                    i++;
                }

                clock_gettime(CLOCK_REALTIME, &t1);
                #ifdef GIN_OMP
                omp_set_num_threads((int)num_threads);
                #endif
                #pragma omp parallel for default(none) shared(graph, i, tasks, num_threads)
                for(int_t k = 0; k < i; k++) {
                    gin_graph_find(graph, tasks[k].str, &tasks[k].matches);
                }
                clock_gettime(CLOCK_REALTIME, &t2);
                utils_find_time += to_sec(t1, t2);

                for(int_t j = 0; j < i; j++) {
                    if(!tasks[j].str) continue;
                    gin_vector_t *decoded_matches = tasks[j].matches;
                    utils_no_matches += decoded_matches->size;
                    for(int_t k = 0; k < decoded_matches->size; k++) {
                        gin_kmer_kv_t *decoded_match = decoded_matches->data[k];
                        fprintf(foutput, "%s:(v:%lld,o:%lld)\n", decoded_match->str->seq, decoded_match->vid, decoded_match->offset);
                    }
                    gin_vector_free(decoded_matches);
                    gin_string_free(tasks[j].str);
                }
                if(exit_flag)
                    break;
            }
            free(buf);
            free(tasks);
            gin_graph_free(graph);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);
            if(verbose) {
                fprintf(stderr, "[gin:utils] Params: Graph file name (-r): %s\n", fref_path);
                if(finput_path) {
                    fprintf(stderr, "[gin:utils] Params: Query file name (-i): %s\n", finput_path);
                }
                fprintf(stderr, "[gin:utils] Params: Read batch size (-b): %lld\n", batch_size);
                fprintf(stderr, "[gin:utils] Params: Threads (-j): %lld\n", num_threads);
                fprintf(stderr, "[gin:utils] Index: Index parse time (s): %lf\n", parse_time);
                fprintf(stderr, "[gin:utils] Decode: Total matches decoded: %lld\n", utils_no_matches);
                fprintf(stderr, "[gin:utils] Decode: Total decoding time (s): %lf\n", utils_find_time);
                fprintf(stderr, "[gin:utils] Decode: Matches decoded per second: %lf\n", (double)utils_no_matches / utils_find_time);
                fprintf(stderr, "[gin:utils] Decode: Time per match decode (s): %lf\n", utils_find_time / (double)utils_no_matches);
                if(utils_queries_processed) {
                    fprintf(stderr, "[gin:utils] Find: Total queries processed: %lld\n",utils_queries_processed);
                    fprintf(stderr, "[gin:utils] Find: Total querying time (s): %lf\n",utils_find_time);
                    fprintf(stderr, "[gin:utils] Find: Queries per second: %lf\n",(double) utils_queries_processed / (double) utils_find_time);
                    fprintf(stderr, "[gin:utils] Find: Time per query (s): %lf\n",(double) utils_find_time / (double) utils_queries_processed);
                }
            }
            break;
        }
        default: {
            fprintf(stderr, "[gin:utils] Unrecognized utils mode, please see gin help utils. Quitting.\n");
            return_code = -1;
            break;
        }
    }
    return return_code;
}

int gin_main_help(gin_mode_t progmode, char *progname) {
    int return_code = 0;
    if(!progname) {
        fprintf(stderr, "%s%s%s\n%s%s",
                "[gin:help] gin! FM-Index like graph indexing algorithm toolkit \n",
                "[gin:help] Needle in a haystack? More like string in a graph. Version 1.1",
                gin_version,
                "[gin:help] Please use gin help <program_name> to learn more about a particular program\n",
                "[gin:help] List of currently available programs: ");
        for(int i = 0; i < gin_mode_no_modes; i++) {
            fprintf(stderr, "%s ", gin_mode_names[i]);
        }
        return 0;
    }

    switch (progmode) {
        case gin_mode_index: {
            fprintf(stderr, "[gin:help] ---------- gin:index ----------\n");
            fprintf(stderr, "[gin:help] gin index indexes a string labelled graph and produces a program specific output for future querying.\n");
            fprintf(stderr, "[gin:help] Input files in rGFA or program specific .ging formats are accepted. Please see ging documentation for further information.\n");
            fprintf(stderr, "[gin:help] Additionally, a vertex permutation can be provided to increase querying times. Please see gin permutation for further information.\n");
            fprintf(stderr, "[gin:help] An identity permutation is used by default.\n");
            fprintf(stderr, "[gin:help] Parameters:\n");
            fprintf(stderr, "\t--input            or -i: Optional parameter. Path to the input file in rGFA or ging format. Default: stdin\n");
            fprintf(stderr, "\t--rgfa             or -g: Optional flag.      Indicates that the input file is an rGFA file. Default: false\n");
            fprintf(stderr, "\t--output           or -o: Optional parameter. Path to the output file, produced in binary gini format. Default: stdout\n");
            fprintf(stderr, "\t--permutation      or -p: Optional parameter. Path to the permutation file. See gin permutation for more help. Default = Identity permutation\n");
            fprintf(stderr, "\t--isa-sample-rate  or -s: Optional parameter. Sampling rate of the suffix array. Reducing this parameter increases query speeds at the cost of larger index files. Default = 32\n");
            fprintf(stderr, "\t--rank-sample-rate or -r: Optional parameter. Frequency of rank caches. Reducing this parameter increases query speeds at the cost of larger index files. Default = 32\n");
            fprintf(stderr, "\t--verbose          or -v: Optional flag.      Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[gin:help] Example invocation: gin index -i mygraph.rgfa -g -o mygraph.gini -p myperm -s 64 -r 64 -v\n");
            return_code = 0;
            break;
        }
        case gin_mode_query : {
            fprintf(stderr, "[gin:help] ---------- gin:query ----------\n");
            fprintf(stderr, "[gin:help] gin query loads a graph index in gini format into memory and runs the provided queries on the graph.\n");
            fprintf(stderr, "[gin:help] Query inputs are expected in the form of one string per line or in FASTQ format.\n");
            fprintf(stderr, "[gin:help] gin query various modes, which are described below:\n");
            fprintf(stderr, "\tcache: Constructs a cache with precomputed suffix array ranges for all string matches up to a specified length. Usage of a cache is highly recommended.\n");
            fprintf(stderr, "\tfind:  Finds matching walk roots. Results are returned as (vertex_id, index) per match per line if --decode is enabled, Otherwise prints the suffix array ranges for each query.\n");
            fprintf(stderr, "[gin:help] Parameters:\n");
            fprintf(stderr, "\t--reference   or -r: Required parameter (find, cache). Path to the index file. See gin index for more help.\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter (find).        Path to the input file containing string queries, with one string per line. Default: stdin\n");
            fprintf(stderr, "\t--fastq       or -f: Optional flag      (find).        Specifies if queries are contained in fastq format. If not in FASTQ format, the program expects on string per line, and the string exit(); to indicate the end of the stream. Default: False\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter (find, cache). Path to the output file. For find, (vertex_id, index) is written to this file if decode is enabled, else suffix array entries are written. For cache, the cache binary is written. Default: stdout\n");
            fprintf(stderr, "\t--cache-depth or -c: Optional parameter (cache).       Specifies the depth of the cache to be constructed. Default: 10\n");
            fprintf(stderr, "\t--cache       or -C: Optional parameter (find).        Path to the index cache. Default: None\n");
            fprintf(stderr, "\t--max-forks   or -m: Optional parameter (find).        Number of maximum forks to be tracked for a query. Setting this to -1 tracks all forks. Default: -1\n");
            fprintf(stderr, "\t--max-matches or -M: Optional parameter (find).        Maximum number of matches decoded for a query. Setting this to -1 decodes all matches. Default: -1\n");
            fprintf(stderr, "\t--decode      or -d: Optional flag      (find).        Decodes the matches into text space. Setting to true may result in combinatorial blowup. Default: False\n");
            fprintf(stderr, "\t--batch-size  or -b: Optional parameter (find).        Number of queries to be read and processed at once. Default: 8\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter (find, cache). Number of threads to be used for parallel querying. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter (find, cache). Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[gin:help] Example invocation (cache): gin query cache -r myindex.gini -o myindex_cache.ginc -j 8 -c 10 -v\n");
            fprintf(stderr, "[gin:help] Example invocation (find):  gin query find  -r myindex.gini -i queries.fastq -f -C myindex_cache.ginc -o results.txt -j 8 -m -1 -M 10 -v\n");

            return_code = 0;
            break;
        }
        case gin_mode_decode : {
            fprintf(stderr, "[gin:help] ---------- gin:decode ----------\n");
            fprintf(stderr, "[gin:help] gin decode loads an bit encoded graph into memory and enumerates full walks from the output of gin query find --decode without -v.\n");
            fprintf(stderr, "[gin:help] Inputs are expected to be in the form of the output of gin query find --decode.\n");
            fprintf(stderr, "[gin:help] gin query has two modes described below:\n");
            fprintf(stderr, "\tencode: Encodes the input ging graph into a program specific format using ceil(log2(s)) bits per character.\n");
            fprintf(stderr, "\twalks:  Enumerates all walks of a given string originating from (vertex,offset) pairs. The output format is of the form:\n"
                            "\t\t<string>:\n"
                            "\t\t\t(o1,oN);v1:...:vN\n"
                            "\t\t\t(o1,oN);v1:...:vN\n");
            fprintf(stderr, "[gin:help] Parameters:\n");
            fprintf(stderr, "\t--reference   or -r: Required parameter (walks). Path to the index file. See gin index for more help.\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter (walks, encode). Path to the input file containing string queries, or the input ging graph to be encoded. Default: stdin\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter (walks, encode). Path to the output file. For encode, the bit encoded graph is written. For walks, resulting walks are written. Default: stdout\n");
            fprintf(stderr, "\t--batch-size  or -b: Optional parameter (walks).         Number of queries to be read and processed at once. Default: 8\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter (walks, encode). Number of threads to be used for parallel querying. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter (walks, encode). Provides more information (time, progress, memory requirements).\n");
            fprintf(stderr, "[gin:help] Example invocation (encode): gin decode encode -i mygraph.ging -o mygraph.gine-v\n");
            fprintf(stderr, "[gin:help] Example invocation (walks):  gin query walks  -r myindex.gine -i queries.query -o results.txt -j 8 -b 16 -v\n");

            return_code = 0;
            break;
        }
        case gin_mode_permutation : {
            fprintf(stderr, "[gin:help] ---------- gin:permutation ----------\n");
            fprintf(stderr, "[gin:help] gin permutation approximates a permutation of vertex indices such that they occur consecutively in the Burrows-Wheeler order.\n");
            fprintf(stderr, "[gin:help] Computing such a permutation is NP-Complete; this program essentially runs simulated annealing and outputs a permutation.\n");
            fprintf(stderr, "[gin:help] Parameters:\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter. Path to the input file in rGFA or ging format. Default: stdin\n");
            fprintf(stderr, "\t--rgfa        or -g: Optional flag.      Indicates that the input file is an rGFA file. Default: false\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter. Path to the output file, one integer as vid per line. Default: stdout\n");
            fprintf(stderr, "\t--permutation or -p: Optional parameter. Path to initial permutation file to start optimizing from. Default: Identity permutation\n");
            fprintf(stderr, "\t--path-span   or -s: Optional flag.      Disallows constraint sets to span multiple vertices. Default: False\n");
            fprintf(stderr, "\t--depth       or -d: Optional parameter. Maximum string length to be considered in the construction of constraint sets. Default: 4\n");
            fprintf(stderr, "\t--temperature or -e: Optional parameter. Sets the initial temperature of the annealing process. Default: 1e6\n");
            fprintf(stderr, "\t--cooling     or -c: Optional parameter. Sets the cooling factor of the annealing process. Default: 0.95\n");
            fprintf(stderr, "\t--time        or -t: Optional parameter. Time after which optimization is terminated in seconds. Default: 15 seconds\n");
            fprintf(stderr, "\t--update      or -u: Optional parameter. Time interval of informative prints in seconds. Default: 3 seconds\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter. Number of threads to be used for parallel cost computation. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter. Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[gin:help] Example invocation: gin permutation -i mygraph.ging -o mygraph-perm.txt -t 300 -u 15 -j 8 -v\n");
            return_code = 0;
            break;
        }
        case gin_mode_utils: {
            fprintf(stderr, "[gin:help] ---------- gin:utils ----------\n");
            fprintf(stderr, "[gin:help] gin utils converts rgfa and FASTQ files to ging files and gin query files, and may also be used to extract the k-mer spectrum of an ging file or find brute force matches.\n");
            fprintf(stderr, "[gin:help] gin utils supports four modes, which are described below:\n");
            fprintf(stderr, "\trgfa2ging:   Converts an rGFA into an ging file. Throws away naming conventions and stable sequences. Everything is assumed to be on the forward strand.\n");
            fprintf(stderr, "\tfastq2query: Extracts the sequence in every read and returns one sequence per line.\n");
            fprintf(stderr, "\tspectrum:    Extracts the k-mer spectrum of an input ging file.\n");
            fprintf(stderr, "\tfind:        Finds (vid,offset) pairs of query matches similar to gin:query, but uses a brute-force approach.\n");
            fprintf(stderr, "[gin:help] Parameters\n");
            fprintf(stderr, "\t--input     or -i: Optional parameter (rgfa2ging, fastq2query, find).          Path to the input file in rGFA (rgfa2ging) or FASTQ (fastq2query) or query (find) format. Default: stdin\n");
            fprintf(stderr, "\t--output    or -o: Optional parameter (rgfa2ging, fastq2query, convert, find). Path to the output file. Default: stdout\n");
            fprintf(stderr, "\t--reference or -r: Requried parameter (find).                                  Path to the input graph file in ging format. Default: stdin\n");
            fprintf(stderr, "\t--threads   or -j: Optional parameter (find).                                  Number of threads.\n") ;
            fprintf(stderr, "\t--verbose   or -v: Optional flag.      Provides more information about the conversion process. Default: false\n");
            fprintf(stderr, "[gin:help] Example invocation: gin utils rgfa2ging -i mygraph.rgfa -o mygraph.ging -v\n");
            fprintf(stderr, "[gin:help] Example invocation: gin utils spectrum -i mygraph.rgfa -o spectrum.txt -v\n");
            fprintf(stderr, "[gin:help] Example invocation: gin utils find -i mygraph.gini -o matches.txt -v -j4 \n");
            return_code = -1;
            break;
        }
        case gin_mode_help: {
            fprintf(stderr, "[gin:help] ---------- gin:help ----------\n");
            fprintf(stderr, "[gin:help] gin help prints general information about how other programs under gin work, then quits.\n");
            fprintf(stderr, "[gin:help] To learn about a particular program, please call gin help <program_name> through the command line.\n");
            fprintf(stderr, "[gin:help] Any other parameters supplied are ignored.\n");
            fprintf(stderr, "[gin:help] List of currently available programs: ");
            for(int i = 0; i < gin_mode_no_modes; i++) {
                fprintf(stderr, "%s ", gin_mode_names[i]);
            }
            fprintf(stderr, "\n");
            return_code = 0;
            break;
        }
        default: {
            fprintf(stderr, "[gin:help] Program %s not recognized. Please run gin help for more information. Quitting.\n",progname);
            return_code = -1;
            break;
        }
    }
    return return_code;
}
