#include "fmd_fmd.h"
#include "fmd_annealing.h"
#include <getopt.h>
#include <time.h>
#include "rgfa_parser.h"
#include "permutation_parser.h"
#include "fmdg_parser.h"
//#include "fmd_fuzzy.h"
#ifdef FMD_OMP
#include <omp.h>
#endif

#define FMD_MAIN_ISA_SAMPLE_RATE_DEFAULT 256
#define FMD_MAIN_RANK_SAMPLE_RATE_DEFAULT 256
#define to_sec(t1,t2) (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9
#define FMD_MAIN_BUF_READ_SIZE 1024

#define FMD_MAIN_QUERY_BUF_LEN 65536
#define FMD_MAIN_QUERY_EXIT_PROMPT "exit();"

char *fmd_version = "1.0";
char* fmd_mode_names[] = {"index","query","permutation","convert","help"};

typedef enum fmd_mode_ {
    fmd_mode_index=0,
    fmd_mode_query=1,
    fmd_mode_permutation=2,
    fmd_mode_convert=3,
    fmd_mode_help=4,
    fmd_mode_no_modes=5
} fmd_mode_t;

char* fmd_query_mode_names[] = {"count", "locate", "enumerate", "breadth"};

typedef enum fmd_query_mode_ {
    fmd_query_mode_count=0,
    fmd_query_mode_locate=1,
    fmd_query_mode_enumerate=2,
    fmd_query_mode_breadth=3,
    fmd_query_mode_no_modes=4,
} fmd_query_mode_t;

char* fmd_convert_mode_names[] = {"rgfa2fmdg", "fastq2query"};

typedef enum fmd_convert_mode_ {
    fmd_convert_mode_rgfa2fmdg=0,
    fmd_convert_mode_fastq2query=1,
    fmd_convert_mode_no_modes=2,
} fmd_convert_mode_t;

int fmd_main_index(int argc, char **argv);
int fmd_main_query(int argc, char **argv, fmd_query_mode_t mode);
int fmd_main_permutation(int argc, char **argv);
int fmd_main_convert(int argc, char **argv, fmd_convert_mode_t mode);
int fmd_main_help(fmd_mode_t progmode, char *progname);

fmd_mode_t fmd_string_to_mode(char *str) {
    if(!str) return fmd_mode_no_modes;
    fmd_mode_t progmode = fmd_mode_no_modes;
    for(int i = 0; i < fmd_mode_no_modes; i++) {
        if(strcmp(str,fmd_mode_names[i]) == 0) {
            progmode = (fmd_mode_t)i;
            break;
        }
    }
    return progmode;
}

fmd_query_mode_t fmd_string_to_qmode(char *str) {
    if(!str) return fmd_query_mode_no_modes;
    fmd_query_mode_t progmode = fmd_query_mode_no_modes;
    for(int i = 0; i < fmd_query_mode_no_modes; i++) {
        if(strcmp(str,fmd_query_mode_names[i]) == 0) {
            progmode = (fmd_query_mode_t)i;
            break;
        }
    }
    return progmode;
}

fmd_convert_mode_t fmd_string_to_cmode(char *str) {
    if(!str) return fmd_convert_mode_no_modes;
    fmd_convert_mode_t progmode = fmd_convert_mode_no_modes;
    for(int i = 0; i < fmd_convert_mode_no_modes; i++) {
        if(strcmp(str,fmd_convert_mode_names[i]) == 0) {
            progmode = (fmd_convert_mode_t)i;
            break;
        }
    }
    return progmode;
}


int main(int argc, char *argv[]) {
    fmd_mode_t progmode = fmd_string_to_mode(argv[1]);
    char **argvp = argv+1;
    int argcp = argc-1;

    int return_code = -1;
    switch (progmode) {
        case fmd_mode_index: {
            return_code = fmd_main_index(argcp, argvp);
            break;
        }
        case fmd_mode_query : {
            char *query_mode_name = argcp > 1 ? argvp[1] : NULL;
            fmd_query_mode_t query_mode = fmd_string_to_qmode(query_mode_name);
            return_code = fmd_main_query(argcp, argvp, query_mode);
            break;
        }
        case fmd_mode_permutation : {
            return_code = fmd_main_permutation(argcp, argvp);
            break;
        }
        case fmd_mode_convert: {
            char *convert_mode_name = argcp > 1 ? argvp[1] : NULL;
            fmd_convert_mode_t convert_mode = fmd_string_to_cmode(convert_mode_name);
            return_code = fmd_main_convert(argcp, argvp,convert_mode);
            break;
        }
        case fmd_mode_help: {
            char *progname = argcp > 1 ? argvp[1] : NULL;
            fmd_mode_t help_mode = fmd_string_to_mode(progname);
            return_code = fmd_main_help(help_mode, progname);
            break;
        }
        default: {
            fprintf(stderr, "[fmd:error] Program %s not recognized. Please run fmd help for more information. Quitting.\n", argv[1]);
            break;
        }
    }
    return return_code;
}

int fmd_main_index(int argc, char **argv) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fperm_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fperm = NULL;
    int_t isa_rate = FMD_MAIN_ISA_SAMPLE_RATE_DEFAULT;
    int_t rank_rate = FMD_MAIN_RANK_SAMPLE_RATE_DEFAULT;
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
                    fprintf(stderr, "[fmd:index] Input path %s could not be opened, quitting.\n", finput_path);
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
                    fprintf(stderr, "[fmd:index] Permutation path %s could not be opened, quitting.\n", fperm_path);
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
                fprintf(stderr, "[fmd:index] Option %s not recognized, please see fmd help index for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    fmd_graph_t *graph = NULL;
    fmd_vector_t *permutation = NULL;
    fmd_fmd_t *fmd = NULL;
    unsigned char *fmd_buf;
    uint64_t fmd_buf_size;

    /**************************************************************************
     * 1 - Parse the input graph the variable graph. By the end of this block,
     * the variable graph should be populated and the files must be closed.
     *************************************************************************/
    if(verbose) {
        fprintf(stderr, "[fmd:index] Parsing input file %s\n", finput_path);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    if(parse_rgfa) {
        rgfa_t *rgfa = rgfa_parse(finput);
        if(finput != stdin) fclose(finput);
        if (!rgfa) {
            fprintf(stderr, "[fmd:index] Failed to parse rGFA file under %s, quitting.\n", finput_path);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);
            if(fperm) fclose(fperm);
            return_code = -1;
            return return_code;
        }
        graph = rgfa_to_fmd_graph(rgfa);
        rgfa_free(rgfa);

    } else {
        graph = fmdg_parse(finput);
        if(!graph) {
            if (finput_path) fclose(finput);
            if (foutput_path) fclose(foutput);
            if (fperm) fclose(fperm);
            fprintf(stderr, "[fmd:index] Malformed fmdg file, quitting.\n");
            return_code = -1;
            return return_code;
        }
    }
    /**************************************************************************
     * 2 - Parse the permutation. By the end of this block, the permutation
     * file must be closed, and the variable permutation should be populated
     *************************************************************************/
    if(verbose) {
        if(fperm_path) fprintf(stderr, "[fmd:index] Input parsed. Parsing permutation file %s\n", fperm_path);
        else fprintf(stderr, "[fmd:index] Input parsed. Using the identity permutation.\n");
    }
    if(fperm) { // can not be provided through stdin; has to be from some file.
        permutation = permutation_parse(fperm);
        if(!permutation) {
            fprintf(stderr, "[fmd:index] Failed to parse permutation file under %s, quitting.\n", fperm_path);
            return_code = -1;
            if(foutput != stdout) fclose(foutput);
            fmd_graph_free(graph);
            return return_code;
        }
        if(permutation->size != graph->vertices->size) {
            fprintf(stderr, "[fmd:index] Permutation cardinality does not match graph vertex size. Quitting.\n");
            return_code = -1;
            if(foutput_path) fclose(foutput);
            fmd_vector_free(permutation);
            fmd_graph_free(graph);
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
        fprintf(stderr, "[fmd:index] Constructing FMD on a graph with:\n");
        fprintf(stderr, "\tNumber of vertices: %lld\n", no_vertices);
        fprintf(stderr, "\tNumber of edges: %lld\n", graph->no_edges);
        int_t no_chars = 0;
        for(int_t i = 0; i < graph->vertex_list->size; i++) {
            void* _;
            fmd_table_lookup(graph->vertices, (void*)i, &_);
            fmd_vertex_t *v = (fmd_vertex_t*)_;
            no_chars += v->label->size;
        }
        fprintf(stderr, "\tAverage in-degree: %.3lf\n", (double)(graph->no_edges) / (double)no_vertices);
        fprintf(stderr, "\tTotal length of vertex labels: %lld\n", no_chars);
        fprintf(stderr, "\tAverage label length per vertex: %.3lf\n", (double)no_chars / (double)no_vertices);
        int_t enc_sz = (no_chars + 2 * no_vertices + fmd_ceil_log2(no_vertices) * no_vertices + 1);
        int_t sa_sz = (int_t)sizeof(int_t) * enc_sz;
        int_t graph_sz = no_chars + ((int_t)sizeof(fmd_vertex_t) + 2*(int_t)sizeof(fmd_vertex_t*)) * no_vertices + (int_t)sizeof(int_t) * 2 * graph->no_edges;
        int_t perm_sz = permutation ? (int_t)sizeof(vid_t) * permutation->size : 0;
        int_t emu = sa_sz + enc_sz + graph_sz + perm_sz;
        fprintf(stderr, "[fmd:index] Estimated peak memory requirement: %lld bytes (%.3lf GB)\n", emu, (double)emu * 1e-9);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    fmd_fmd_init(&fmd,
                 graph,
                 permutation,
                 FMD_FMD_DEFAULT_c_0,
                 FMD_FMD_DEFAULT_c_1,
                 rank_rate,
                 isa_rate);
    clock_gettime(CLOCK_REALTIME, &t2);
    index_time = to_sec(t1,t2);
    fmd_graph_free(graph);
    fmd_vector_free(permutation);
    /**************************************************************************
    * 4 - Convert the FM-Index to a buffer and write to disk
    *************************************************************************/
    clock_gettime(CLOCK_REALTIME, &t1);
    fmd_fmd_serialize_to_buffer(fmd, &fmd_buf, &fmd_buf_size);
    fmd_fmd_free(fmd);
    if(!fmd_buf || !fmd_buf_size) {
        fprintf(stderr, "[fmd:index] Could not convert graph into a bitstream, please send a PR. Quitting.\n");
        return_code = -1;
        if(foutput != stdout) fclose(foutput);
        return return_code;
    }

    if(foutput_path) {
        foutput = fopen(foutput_path, "w");
        if(!foutput) {
            fprintf(stderr, "[fmd:index] Output path %s could not be opened, quitting.\n", foutput_path);
            free(fmd_buf);
            return_code = -1;
            return return_code;
        }
    }

    fwrite(fmd_buf,sizeof(unsigned char),fmd_buf_size,foutput);
    if(foutput_path) fclose(foutput);
    clock_gettime(CLOCK_REALTIME, &t2);
    write_time = to_sec(t1,t2);
    if(verbose) {
        fprintf(stderr, "[fmd:index] Resulting index size: %lld bytes (%.3lf GB).\n", fmd_buf_size, (double) fmd_buf_size * 1e-9);
        fprintf(stderr, "[fmd:index] Timings in seconds: \n");
        fprintf(stderr, "\t Parsing:  %lf\n", parse_time);
        fprintf(stderr, "\t Indexing: %lf\n", index_time);
        fprintf(stderr, "\t Writing:  %lf\n", write_time);
        fprintf(stderr, "[fmd:index] Total relevant runtime: %lf seconds.\n", parse_time + index_time + write_time);
    }

    free(fmd_buf);
    return return_code;
}

int fmd_main_query(int argc, char **argv, fmd_query_mode_t mode) {
    int return_code = 0;
    if(mode == fmd_query_mode_no_modes) {
        fprintf(stderr, "[fmd:query] Invalid query mode %s supplied. Please run fmd help query for more information. Quitting...\n", argv[0]);
        return_code = -1;
        return return_code;
    }

    char *fref_path = NULL;
    char *finput_path = NULL;
    char *foutput_path = NULL;

    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fref = NULL;

    int_t num_threads = 1;
    int_t batch_size = 8;
    int_t max_matches = -1;
    bool parse_fastq = false;
    bool verbose = false;

    struct timespec t1;
    struct timespec t2;

    // statistics for path enumeration
    int_t queries_processed = 0;
    double query_time = .0;

    int_t no_matching_forks = 0;
    int_t no_missing_forks = 0;
    int_t no_matching_count = 0; // to compare with DFS
    int_t no_multiple_vertex_span_forks = 0; // keep track of vertices spanning multiple vertices
    int_t no_multiple_vertex_span_matches = 0; // same but keeps the count

    static struct option options[] = {
            {"reference",  required_argument, NULL, 'r'},
            {"input",      required_argument, NULL, 'i'},
            {"fastq",      no_argument,       NULL, 'f'},
            {"output",     required_argument, NULL, 'o'},
            {"max-matches",required_argument, NULL, 'm'},
            {"batch-size", required_argument, NULL, 'b'},
            {"threads",    required_argument, NULL, 'j'},
            {"verbose",    no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "r:i:fo:m:b:j:v", options, &optindex)) != -1) {
        switch(c) {
            case 'r': {
                fref_path = optarg;
                break;
            }
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[fmd:query] Input query path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'f': {
                parse_fastq = true;
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'm': {
                max_matches = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'b': {
                batch_size = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'j': {
                num_threads = (int_t)strtoull(optarg, NULL, 10);
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[fmd:query] Option %s not recognized, please see fmd help query for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    /**************************************************************************
    * 1 - Read the graph index from disk, fmd_buf and fmd_buf_size should be
    * populated and fref be closed by the end of this block
    **************************************************************************/

    fref = fopen(fref_path, "r");
    if(!fref) {
        fprintf(stderr, "[fmd:query] Could not open index reference file %s, quitting.\n", fref_path);
        return_code = -1;
        return return_code;
    }

    uint64_t fmd_buf_capacity = 65536;
    uint64_t fmd_buf_size = 0;
    uint8_t *fmd_buf = calloc(fmd_buf_capacity, sizeof(uint8_t));

    uint64_t read_size;
    uint8_t *read_buf = calloc(FMD_MAIN_BUF_READ_SIZE, sizeof(uint8_t));
    while((read_size = fread(read_buf,sizeof(uint8_t), FMD_MAIN_BUF_READ_SIZE, fref))) {
        if(read_size + fmd_buf_size > fmd_buf_capacity) {
            fmd_buf = realloc(fmd_buf, (fmd_buf_capacity *= 2)*sizeof(uint8_t));
            memset(fmd_buf + (fmd_buf_capacity/2), 0, sizeof(uint8_t)*(fmd_buf_capacity/2));
        }
        memcpy(fmd_buf+fmd_buf_size,read_buf,read_size);
        fmd_buf_size+=read_size;
    }
    fmd_buf = realloc(fmd_buf, fmd_buf_size*sizeof(uint8_t));
    free(read_buf);
    fclose(fref);

    /**************************************************************************
    * 2 - Parse the bitstream into data structures, by the end of this block
    * fmd should be populated and fmd_buf should be freed
    **************************************************************************/
    if(verbose) {
        fprintf(stderr, "[fmd:query] Parsing reference index into data structures.\n");
    }
    fmd_fmd_t *fmd;
    fmd_fmd_serialize_from_buffer(&fmd, fmd_buf, fmd_buf_size);
    if(!fmd) {
        fprintf(stderr, "[fmd:query] Error encountered while parsing reference index. Quitting.\n");
        return_code = -1;
        return return_code;
    }
    /**************************************************************************
    * 3 - Parse queries from the input stream and query depending on the mode
    **************************************************************************/
    if(verbose) {
        fprintf(stderr, "[fmd:query] Parsing and launching queries with %d threads and %d batch size.\n", (int)num_threads, (int)batch_size);
    }
    if(finput_path) {
        finput = fopen(finput_path, "r");
        if(!finput) {
            fprintf(stderr,"[fmd:query] Can't open input file %s, quitting.\n", finput_path);
            return_code = -1;
            return return_code;
        }
    }
    if(parse_fastq) {
        // not yet implemented
        fprintf(stderr, "[fmd:query] Fastq parsing mode not yet implemented. Quitting.\n");
        return_code = -1;
        fmd_fmd_free(fmd);
        return return_code;
    } else {
#ifdef FMD_OMP
        int_t i;
        char *buf = calloc(FMD_MAIN_QUERY_BUF_LEN, sizeof(char*));
        typedef struct query_task_ {
            fmd_string_t *str;
            uint64_t count;
            fmd_vector_t *paths_or_locs;
            fmd_vector_t *partial_matches;
        } query_task_t;
        omp_set_num_threads((int)num_threads);
        query_task_t *tasks = calloc(batch_size, sizeof(query_task_t));
        bool exit_flag = false;

        while(true) {
            i = 0;
            while (i < batch_size && fgets(buf, FMD_MAIN_QUERY_BUF_LEN, finput)) {
                int_t len = (int_t)strlen(buf);
                if(buf[len-1] == '\n') buf[len-1] = 0; // get rid of the end line character
                if(!strcmp(buf,FMD_MAIN_QUERY_EXIT_PROMPT)) {
                    exit_flag = true;
                    break;
                }
                fmd_string_init_cstr(&tasks[i].str, buf);
                queries_processed++;
                i++;
            }

            clock_gettime(CLOCK_REALTIME, &t1);
            omp_set_num_threads((int)num_threads);
            switch(mode) {
                case fmd_query_mode_count: {
                    #pragma omp parallel for default(none) shared(fmd, i, tasks)
                    for(int_t k = 0; k < i; k++) {
                        tasks[k].count = fmd_fmd_query_count(fmd, tasks[k].str);
                    }
                    break;
                }
                case fmd_query_mode_locate: {
                    #pragma omp parallel for default(none) shared(fmd, i, tasks)
                    for(int_t k = 0; k < i; k++) {
                        tasks[k].paths_or_locs = fmd_fmd_query_locate_basic(fmd, tasks[k].str);
                    }
                    break;
                }
                case fmd_query_mode_enumerate: {
                    #pragma omp parallel for default(none) shared(fmd, i, tasks, num_threads, max_matches)
                    for(int_t k = 0; k < i; k++) {
                        fmd_fmd_query_locate_paths_omp(fmd, tasks[k].str, max_matches, &tasks[k].paths_or_locs, &tasks[k].partial_matches, num_threads);
                    }
                    break;
                }
                case fmd_query_mode_breadth: {
                    #pragma omp parallel for default(none) shared(fmd, i, tasks, num_threads, max_matches)
                    for(int_t k = 0; k < i; k++) {
                        fmd_fmd_locate_paths_breadth_first(fmd, tasks[k].str, max_matches, &tasks[k].paths_or_locs, &tasks[k].partial_matches);
                    }
                    break;
                }
                default: {
                    break;
                }
            }
            clock_gettime(CLOCK_REALTIME, &t2);
            query_time += to_sec(t1,t2);

            switch (mode) {
                case fmd_query_mode_count: {
                    for (int_t j = 0; j < i; j++) {
                        if (!tasks[j].str) continue;
                        fprintf(foutput, "%lld\n", tasks[j].count);
                        fmd_string_free(tasks[j].str);
                        tasks[j].str = NULL;
                    }
                    break;
                }
                case fmd_query_mode_locate: {
                    for (int_t j = 0; j < i; j++) {
                        if (!tasks[j].str) continue;
                        for (int_t k = 0; k < tasks[j].paths_or_locs->size - 1; k++) {
                            fprintf(foutput, "%lld, ", (uint64_t) tasks[j].paths_or_locs->data[k]);
                        }
                        if (tasks[j].paths_or_locs->size) fprintf(foutput, "%lld, ",
                                                                (uint64_t) tasks[j].paths_or_locs->data[
                                                                        tasks[j].paths_or_locs->size - 1]);
                        else fprintf(foutput, "-\n");
                        fmd_string_free(tasks[j].str);
                        fmd_vector_free(tasks[j].paths_or_locs);
                        tasks[j].str = NULL;
                    }
                    break;
                }
                case fmd_query_mode_breadth: {
                    for (int_t j = 0; j < i; j++) {
                        if (!tasks[j].str) continue;
                        fmd_vector_t *match_lists;
                        fmd_fmd_query_locate_paths_topologise(&match_lists, tasks[j].str, tasks[j].paths_or_locs);
                        no_matching_forks += tasks[j].paths_or_locs->size;
                        no_missing_forks += tasks[j].partial_matches->size;
                        for (int_t k = 0; k < match_lists->size; k++) {
                            fmd_fmd_match_list_t *list = match_lists->data[k];
                            fmd_fmd_match_node_t *root = (fmd_fmd_match_node_t *) list->head;
                            no_matching_count += root->sa_hi - root->sa_lo;
                            no_multiple_vertex_span_forks += (root->v_lo > -1);
                            no_multiple_vertex_span_matches += (root->v_hi > root->v_lo) * (root->sa_hi - root->sa_lo);
                            fprintf(foutput, "(%s,v:(%lld,%lld),sa:(%lld,%lld))", root->matching_substring->seq, root->v_lo,
                                    root->v_hi, root->sa_lo, root->sa_hi);
                            root = (fmd_fmd_match_node_t *) root->next;
                            while (root != list->dummy) {
                                fprintf(foutput, "->(%s,v:(%lld,%lld),sa:(%lld,%lld))", root->matching_substring->seq, root->v_lo,
                                        root->v_hi, root->sa_lo, root->sa_hi);
                                root = (fmd_fmd_match_node_t *) root->next;
                            }
                            if (!verbose) fprintf(foutput, "\n");
                            else fprintf(foutput, ": %s\n", tasks[j].str->seq);
                        }
                        if (!tasks[j].paths_or_locs->size) fprintf(foutput, "-\n");
                        fprintf(foutput, "\n");
                        fmd_fmd_query_locate_paths_topologise_free(match_lists);
                        fmd_fmd_locate_paths_result_free(tasks[j].paths_or_locs, tasks[j].partial_matches);
                        fmd_string_free(tasks[j].str);
                        tasks[j].str = NULL;
                    }
                    break;
                }
                case fmd_query_mode_enumerate:{
                    for (int_t j = 0; j < i; j++) {
                        if (!tasks[j].str) continue;
                        no_matching_forks += tasks[j].paths_or_locs->size;
                        no_missing_forks += tasks[j].partial_matches->size;
                        for (int_t k = 0; k < tasks[j].paths_or_locs->size; k++) {
                            fmd_fork_node_t *root = (fmd_fork_node_t *) tasks[j].paths_or_locs->data[k];
                            no_matching_count += root->sa_hi - root->sa_lo;
                            no_multiple_vertex_span_forks += (root->vertex_lo > -1);
                            no_multiple_vertex_span_matches += (root->vertex_hi > root->vertex_lo) * (root->sa_hi - root->sa_lo);
                            fprintf(foutput, "%s:(v:(%lld,%lld),sa:(%lld,%lld),pos:%lld)", tasks[j].str->seq, root->vertex_lo,
                                    root->vertex_hi, root->sa_lo, root->sa_hi, root->pos);
                            root = (fmd_fork_node_t *) root->parent;
                            while (root) {
                                fprintf(foutput, "->(v:(%lld,%lld),sa:(%lld,%lld),pos:%lld)", root->vertex_lo,
                                        root->vertex_hi, root->sa_lo, root->sa_hi, root->pos);
                                root = (fmd_fork_node_t *) root->parent;
                            }
                            if (!verbose) fprintf(foutput, "\n");
                            else fprintf(foutput, ": %s\n", tasks[j].str->seq);
                        }
                        if (!tasks[j].paths_or_locs->size) fprintf(foutput, "-\n");
                        fprintf(foutput, "\n");
                        fmd_fmd_locate_paths_result_free(tasks[j].paths_or_locs, tasks[j].partial_matches);
                        fmd_string_free(tasks[j].str);
                        tasks[j].str = NULL;
                    }
                    break;
                }
                default: {
                    break;
                }
            }
            if(exit_flag)
                break;
        }
        free(buf);
        free(tasks);
#else
        char *query_cstr = calloc(FMD_MAIN_QUERY_BUF_LEN, sizeof(char));
        clock_gettime(CLOCK_REALTIME, &t1);
        while(fgets(query_cstr, FMD_MAIN_QUERY_BUF_LEN, finput)) {
            size_t query_len = strlen(query_cstr);
            if(!query_len) continue; // skip empty lines
            if(query_cstr[query_len-1] == '\n')query_cstr[query_len-1] = 0; // get rid of the line break
            if(!strcmp(query_cstr,FMD_MAIN_QUERY_EXIT_PROMPT)) break; // break if exit prompt is given

            fmd_string_t *query;
            fmd_string_init_cstr(&query, query_cstr);
            queries_processed++;

            switch(mode) {
                case fmd_query_mode_count: {
                    uint64_t count = fmd_fmd_query_count(fmd, query);
                    fprintf(foutput,"%lld\n",count);
                    break;
                }
                case fmd_query_mode_locate: {
                    fmd_vector_t *locs = fmd_fmd_query_locate_basic(fmd, query);
                    for(int_t i = 0; i < locs->size-1; i++) {
                        fprintf(foutput, "%lld, ", (uint64_t)locs->data[i]);
                    }
                    if(locs->size) fprintf(foutput, "%lld, ", (uint64_t)locs->data[locs->size-1]);
                    free(locs);
                    break;
                }
                case fmd_query_mode_enumerate: {
                    fmd_vector_t *paths, *graveyard;
                    fmd_fmd_query_locate_paths(fmd, query, max_matches, &paths, &graveyard);
                    no_matching_forks += paths->size;
                    no_missing_forks += graveyard->size;
                    for(int_t i = 0; i < paths->size; i++) {
                        fmd_fork_node_t *root = (fmd_fork_node_t*)paths->data[i];
                        fprintf(foutput, "(v:(%lld,%lld),sa:(%lld,%lld),pos:%lld)", root->vertex_lo, root->vertex_hi, root->sa_lo, root->sa_hi, root->pos);
                        root = (fmd_fork_node_t*)root->parent;
                        no_matching_count += root->sa_hi - root->sa_lo;
                        no_multiple_vertex_span_forks += (root->vertex_lo > -1);
                        no_multiple_vertex_span_matches += (root->vertex_hi > root->vertex_lo) * (root->sa_hi - root->sa_lo);
                        while(root) {
                            fprintf(foutput, "->(v:(%lld,%lld),sa:(%lld,%lld),pos:%lld)", root->vertex_lo, root->vertex_hi, root->sa_lo, root->sa_hi, root->pos);
                            root = (fmd_fork_node_t*)root->parent;
                        }
                        fprintf(foutput, "\n");
                    }
                    fprintf(foutput, "\n");
                    fmd_fmd_locate_paths_result_free(paths, graveyard);
                }
                default: {
                    break;
                }
            }
            fmd_string_free(query);
        }
        clock_gettime(CLOCK_REALTIME, &t2);
        query_time = to_sec(t1,t2);
        free(query_cstr);
#endif
    }
    /**************************************************************************
    * 4 - Free all data structures and return
    **************************************************************************/
    if(finput_path) fclose(finput);
    if(foutput_path) fclose(foutput);

    if(verbose) {
        fprintf(stderr, "[fmd:query] Total querying time in seconds: %lf\n",query_time);
        fprintf(stderr, "[fmd:query] Number of queries processed: %lld\n",queries_processed);
        if(queries_processed)
            fprintf(stderr, "[fmd:query] Average time per query: %lf\n",(double)query_time / (double)queries_processed);
        if(mode == fmd_query_mode_enumerate && queries_processed && no_matching_forks) {
            fprintf(stderr, "[fmd:query] Forks:\n");
            fprintf(stderr, "[fmd:query] Number of matching forks: %lld\n",no_matching_forks);
            fprintf(stderr, "[fmd:query] Number of partial forks: %lld\n",no_missing_forks);
            fprintf(stderr, "[fmd:query] Number of forks spanning multiple vertices: %lld\n", no_multiple_vertex_span_forks);
            fprintf(stderr, "[fmd:query] Matches:\n");
            fprintf(stderr, "[fmd:query] Number of matches: %lld\n",no_matching_count);
            fprintf(stderr, "[fmd:query] Number of matches spanning multiple vertices: %lld\n", no_multiple_vertex_span_matches);
            fprintf(stderr, "[fmd:query] Aggregate statistics:\n");
            fprintf(stderr, "[fmd:query] Average matches per matching fork: %.6lf\n",(double)no_matching_count / (double)no_matching_forks);
            fprintf(stderr, "[fmd:query] Average matches per fork: %.6lf\n",(double)no_matching_count / ((double)no_matching_forks + (double)no_missing_forks));
            fprintf(stderr, "[fmd:query] Average forks per query: %.6lf\n",((double)no_matching_forks + (double)no_missing_forks) / (double)queries_processed);
            fprintf(stderr, "[fmd:query] Average matching forks per query: %.6lf\n",((double)no_matching_forks) / (double)queries_processed);
            fprintf(stderr, "[fmd:query] Aggregate timings:\n");
            fprintf(stderr, "[fmd:query] Average time per fork: %.6lf\n", (double)query_time / ((double)no_matching_forks + (double)no_missing_forks));
            fprintf(stderr, "[fmd:query] Average time per matching fork: %.6lf\n", (double)query_time / ((double)no_matching_forks));
            fprintf(stderr, "[fmd:query] Average time per match: %.8lf\n", (double)query_time / ((double)no_matching_count));
        }
    }

    fmd_fmd_free(fmd);
    return return_code;
}

int fmd_main_permutation(int argc, char **argv) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    char *fperm_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;
    FILE *fperm = NULL;

    bool parse_rgfa = false;
    bool multiple_vertex_span = false;

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
                    fprintf(stderr, "[fmd:permutation] Input path %s could not be opened, quitting.\n", finput_path);
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
                    fprintf(stderr, "[fmd:permutation] Permutation path %s could not be opened, quitting.\n", fperm_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 's' : {
                multiple_vertex_span = true;
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
                fprintf(stderr, "[fmd:permutation] Option %s not recognized, please see fmd help permutation for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    fmd_graph_t *graph = NULL;
    fmd_vector_t *permutation = NULL;
    fmd_vector_t *constraints = NULL;
    fmd_annealing_t *ann = NULL;
    fmd_vector_t *optimized_permutation = NULL;

    /**************************************************************************
     * 1 - Parse the input graph the variable graph. By the end of this block,
     * the variable graph should be populated and the files must be closed.
     *************************************************************************/
    if(verbose) {
        fprintf(stderr, "[fmd:permutation] Parsing input file %s\n", finput_path);
    }
    clock_gettime(CLOCK_REALTIME, &t1);
    if(parse_rgfa) {
        rgfa_t *rgfa = rgfa_parse(finput);
        if(finput != stdin) fclose(finput);
        if (!rgfa) {
            fprintf(stderr, "[fmd:permutation] Failed to parse rGFA file under %s, quitting.\n", finput_path);
            if(finput_path) fclose(finput);
            if(foutput_path) fclose(foutput);
            if(fperm) fclose(fperm);
            return_code = -1;
            return return_code;
        }
        graph = rgfa_to_fmd_graph(rgfa);
        rgfa_free(rgfa);

    } else {
        graph = fmdg_parse(finput);
        if(!graph) {
            if (finput_path) fclose(finput);
            if (foutput_path) fclose(foutput);
            if (fperm) fclose(fperm);
            fprintf(stderr, "[fmd:permutation] Malformed fmdg file, quitting.\n");
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
            fprintf(stderr, "[fmd:permutation] Failed to parse permutation file under %s, quitting.\n", fperm_path);
            return_code = -1;
            if(foutput != stdout) fclose(foutput);
            fmd_graph_free(graph);
            return return_code;
        }
        if(permutation->size != graph->vertices->size) {
            fprintf(stderr, "[fmd:permutation] Permutation cardinality does not match graph vertex size. Quitting.\n");
            return_code = -1;
            if(foutput_path) fclose(foutput);
            fmd_vector_free(permutation);
            fmd_graph_free(graph);
            return return_code;
        }
    }

    /**************************************************************************
     * 3 - Enumerate constraint sets, by the end of this block constraints
     * must be populated with constraint sets.
     **************************************************************************/
    if(verbose) {
        fprintf(stderr, "[fmd:permutation] Extracting constraint sets for each occurring prefix.\n");
        char *spanstr = multiple_vertex_span ? "True" : "False";
        fprintf(stderr, "[fmd:permutation] Constraint sets can span multiple vertices: %s\n", spanstr);
    }
    fmd_constraint_set_enumerate(&constraints, graph, depth, multiple_vertex_span);
    if(!constraints) {
        fprintf(stderr, "[fmd:permutation] Something went wrong during constraint set extraction. Quitting.\n");
        fmd_graph_free(graph);
        return_code = -1;
        return return_code;
    }

    /**************************************************************************
     * 4 - Configure the annealing optimizer and begin optimizing
     *************************************************************************/
#ifdef FMD_OMP
    omp_set_num_threads((int)num_threads);
#endif
    fmd_annealing_configure(&ann, graph, constraints, permutation, temperature, 1, cooling_factor, 0);
    fmd_graph_free(graph);
    fmd_vector_free(constraints);
    if(verbose) {
        fprintf(stderr, "[fmd:permutation] Optimization begins with initial cost %.3lf for depth=%lld for %lld seconds.\n", ann->cur_cost, depth, time);
    }
    int_t no_iterations = 1 + (time - 1) / update;
    for(int_t i = 0; i < no_iterations; i++) {
        fmd_annealing_iterate_seconds(ann, (int_t)update);
        if(verbose) {
            fprintf(stderr, "\tIteration %d, best_cost = %.3lf, cur_cost = %.3lf\n", ann->cur_iter + 1, ann->best_cost_so_far, ann->cur_cost);
        }
    }

    /**************************************************************************
     * 5 - Cleanup and write to the output
     *************************************************************************/
    fmd_annealing_get_permutation(ann, &optimized_permutation);
    //fmd_annealing_free(ann);
    if(foutput_path) {
        foutput = fopen(foutput_path, "w");
        if(!foutput) {
            fprintf(stderr, "[fmd:index] Output path %s could not be opened, quitting.\n", foutput_path);
            return_code = -1;
            return return_code;
        }
    }

    for(int_t i = 0; i < optimized_permutation->size; i++) {
        fprintf(foutput, "%lld\n", (int_t)optimized_permutation->data[i]);
    }

    if(foutput_path) fclose(foutput);
    fmd_vector_free(permutation);
    fmd_vector_free(optimized_permutation);

    return return_code;
}

int fmd_main_convert(int argc, char **argv, fmd_convert_mode_t mode) {
    char *finput_path = NULL;
    char *foutput_path = NULL;
    FILE *finput = stdin;
    FILE *foutput = stdout;

    bool verbose = false;

    struct timespec t1;
    struct timespec t2;

    double parse_time = 0.0;
    double convert_time = 0.0;
    double write_time = 0.0;

    int return_code = 0;
    static struct option options[] = {
            {"input",       required_argument, NULL, 'i'},
            {"output",      required_argument, NULL, 'o'},
            {"verbose",     no_argument,       NULL, 'v'},
    };
    opterr = 0;
    int optindex,c;
    while((c = getopt_long(argc, argv, "i:o:v", options, &optindex)) != -1) {
        switch(c) {
            case 'i': {
                finput_path = optarg;
                finput = fopen(finput_path, "r");
                if(!finput) {
                    fprintf(stderr, "[fmd:convert] Input path %s could not be opened, quitting.\n", finput_path);
                    return_code = -1;
                    return return_code;
                }
                break;
            }
            case 'o': {
                foutput_path = optarg;
                break;
            }
            case 'v': {
                verbose = true;
                break;
            }
            default: {
                fprintf(stderr, "[fmd:convert] Option %s not recognized, please see fmd help convert for more.\n",optarg);
                return_code = -1;
                break;
            }
        }
    }
    fmd_graph_t *graph = NULL;
    switch (mode) {
        case fmd_convert_mode_rgfa2fmdg: {
            // first parse the graph
            clock_gettime(CLOCK_REALTIME, &t1);
            if(verbose) {
                if(finput_path) fprintf(stderr, "[fmd:convert] Parsing input file %s\n", finput_path);
                else fprintf(stderr, "[fmd:convert] Parsing input from stdin.\n");
            }
            rgfa_t *rgfa = rgfa_parse(finput);
            clock_gettime(CLOCK_REALTIME, &t2);
            parse_time = to_sec(t1,t2);
            if(finput_path) fclose(finput);
            if(!rgfa) {
                fprintf(stderr, "[fmd:convert] Failed to parse rgfa file quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(verbose) {
                if(finput_path) fprintf(stderr, "[fmd:convert] Converting to fmdg format.\n");
            }
            clock_gettime(CLOCK_REALTIME, &t1);
            graph = rgfa_to_fmd_graph(rgfa);
            clock_gettime(CLOCK_REALTIME, &t2);
            convert_time = to_sec(t1,t2);
            rgfa_free(rgfa);
            if(!graph) {
                fprintf(stderr, "[fmd:convert] Failed to convert rgfa_t to fmd_graph_t. Quitting.\n");
                return_code = -1;
                return return_code;
            }
            if(foutput_path) {
                foutput = fopen(foutput_path, "w");
                if(!foutput) {
                    fprintf(stderr, "[fmd:convert] Failed to open output path %s, quitting\n", foutput_path);
                    fmd_graph_free(graph);
                    return_code = -1;
                    return return_code;
                }
            }
            if(verbose) {
                if(foutput_path) fprintf(stderr, "[fmd:convert] Writing converted graph to %s.\n", foutput_path);
                else fprintf(stderr, "[fmd:convert] Writing converted graph to stdout.\n");
            }
            clock_gettime(CLOCK_REALTIME, &t1);
            fmdg_write(foutput, graph);
            clock_gettime(CLOCK_REALTIME, &t2);
            write_time = to_sec(t1,t2);
            if(foutput_path) fclose(foutput);
            if(verbose) {
                fprintf(stderr, "[fmd:convert] Timings in seconds: \n");
                fprintf(stderr, "\t Parsing:    %.3lf\n", parse_time);
                fprintf(stderr, "\t Converting: %.3lf\n", convert_time);
                fprintf(stderr, "\t Writing:    %.3lf\n", write_time);
                fprintf(stderr, "[fmd:index] Total relevant runtime: %lf seconds.\n", parse_time + convert_time + write_time);
            }
            break;
        }
        case fmd_convert_mode_fastq2query: {
            fprintf(stderr, "[fmd:convert] fasta2query not implemented yet. Quitting.\n");
            return_code = -1;
            break;
        }
        default: {
            fprintf(stderr, "[fmd:convert] Unrecognized conversion mode, please see fmd help convert. Quitting.\n");
            return_code = -1;
            break;
        }
    }
    return return_code;
}

int fmd_main_help(fmd_mode_t progmode, char *progname) {
    int return_code = 0;
    if(!progname) {
        fprintf(stderr, "%s%s%s\n%s%s",
                "[fmd:help] fmd! FM-Index like graph indexing algorithm toolkit \n",
                "[fmd:help] Needle in a haystack? More like string in a graph. Version 1.0",
                fmd_version,
                "[fmd:help] Please use fmd help <program_name> to learn more about a particular program\n",
                "[fmd:help] List of currently available programs: ");
        for(int i = 0; i < fmd_mode_no_modes; i++) {
            fprintf(stderr, "%s ", fmd_mode_names[i]);
        }
        return 0;
    }

    switch (progmode) {
        case fmd_mode_index: {
            fprintf(stderr, "[fmd:help] ---------- fmd:index ----------\n");
            fprintf(stderr, "[fmd:help] fmd index indexes a string labelled graph and produces a program specific output for future querying.\n");
            fprintf(stderr, "[fmd:help] Input files in rGFA or program specific .fmdg formats are accepted. Please see fmdg documentation for further information.\n");
            fprintf(stderr, "[fmd:help] Additionally, a vertex permutation can be provided to increase querying times. Please see fmd permutation for further information.\n");
            fprintf(stderr, "[fmd:help] An identity permutation is used by default.\n");
            fprintf(stderr, "[fmd:help] Parameters:\n");
            fprintf(stderr, "\t--input            or -i: Optional parameter. Path to the input file in rGFA or fmdg format. Default: stdin\n");
            fprintf(stderr, "\t--rgfa             or -g: Optional flag. Indicates that the input file is an rGFA file. Default: false\n");
            fprintf(stderr, "\t--output           or -o: Optional parameter. Path to the output file, produced in binary fmdi format. Default: stdout\n");
            fprintf(stderr, "\t--permutation      or -p: Optional parameter. Path to the permutation file. See fmd permutation for more help. Default = Identity permutation\n");
            fprintf(stderr, "\t--isa-sample-rate  or -s: Optional parameter. Sampling rate of the suffix array. Reducing this parameter increases query speeds at the cost of larger index files. Default = 256\n");
            fprintf(stderr, "\t--rank-sample-rate or -r: Optional parameter. Frequency of rank caches. Reducing this parameter increases query speeds at the cost of larger index files. Default = 256\n");
            fprintf(stderr, "\t--verbose          or -v: Optional flag.      Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd index -i mygraph.rgfa -g -o mygraph.fmdi -p myperm -s 64 -r 64 -v\n");
            return_code = 0;
            break;
        }
        case fmd_mode_query : {
            fprintf(stderr, "[fmd:help] ---------- fmd:query ----------\n");
            fprintf(stderr, "[fmd:help] fmd query loads a graph index in fmdi format into memory and runs the provided queries on the graph.\n");
            fprintf(stderr, "[fmd:help] Query inputs are expected in the form of one string per line or in FASTQ format.\n");
            fprintf(stderr, "[fmd:help] fmd query supports many querying modes, which are described below:\n");
            fprintf(stderr, "\tcount:     Counts the occurrences of the string in the space of strings induced by the string graph. Results are returned as a single integer per line.\n");
            fprintf(stderr, "\tlocate:    Locates the vertex IDs and positions into vertices of the string matches. Results are returned as (vertex_id, index) tuples per line.\n");
            fprintf(stderr, "\tenumerate: Enumerates the paths to which vertices match. Results are returned as (vertex_id, index) (vertex_id) ... (vertex_id) (vertex_id index) per match per line.\n");
            fprintf(stderr, "[fmd:help] Parameters:\n");
            fprintf(stderr, "\t--reference   or -r: Required parameter. Path to the index file. See fmd index for more help.\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter. Path to the input file containing string queries, with one string per line. Default: stdin\n");
            fprintf(stderr, "\t--fastq       or -f: Optional flag.      Specifies if queries are contained in fastq format. Default: False\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter. Path to the output file, produced in one of the query mode formats described above. Default: stdout\n");
            fprintf(stderr, "\t--max-matches or -m: Optimal parameter.  Number of maximum matches to be returned for a query. Setting this to -1 returns all matches. Default: -1\n");
            fprintf(stderr, "\t--batch-size  or -b: Optional parameter. Number of queries to be read and processed at once. Default: 8\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter. Number of threads to be used for parallel querying. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter. Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd query enumerate -r myindex.fmdi -i queries.fastq -f -o results.txt -j 8 -m 10 -v\n");
            return_code = 0;
            break;
        }
        case fmd_mode_permutation : {
            fprintf(stderr, "[fmd:help] ---------- fmd:permutation ----------\n");
            fprintf(stderr, "[fmd:help] fmd permutation approximates a permutation of vertex indices such that they occur consecutively in the Burrows-Wheeler order.\n");
            fprintf(stderr, "[fmd:help] Computing such a permutation is NP-Complete; this program essentially runs simulated annealing and outputs a permutation.\n");
            fprintf(stderr, "[fmd:help] Parameters:\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter. Path to the input file in rGFA or fmdg format. Default: stdin\n");
            fprintf(stderr, "\t--rgfa        or -g: Optional flag. Indicates that the input file is an rGFA file. Default: false\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter. Path to the output file, one integer as vid per line. Default: stdout\n");
            fprintf(stderr, "\t--permutation or -p: Optional parameter. Path to initial permutation file to start optimizing from. Default: Identity permutation\n");
            fprintf(stderr, "\t--path-span   or -s: Optional flag. Allows constraint sets to span multiple vertices. Default: False\n");
            fprintf(stderr, "\t--depth       or -d: Optional parameter. Maximum string length to be considered in the construction of constraint sets. Default: 4\n");
            fprintf(stderr, "\t--temperature or -e: Optional parameter. Sets the initial temperature of the annealing process. Default: 1e6\n");
            fprintf(stderr, "\t--cooling     or -c: Optional parameter. Sets the cooling factor of the annealing process. Default: 0.95\n");
            fprintf(stderr, "\t--time        or -t: Optional parameter. Time after which optimization is terminated in seconds. Default: 15 seconds\n");
            fprintf(stderr, "\t--update      or -u: Optional parameter. Time interval of informative prints in seconds. Default: 3 seconds\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter. Number of threads to be used for parallel cost computation. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter. Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd permutation -i mygraph.fmdg -o mygraph-perm.txt -t 300 -u 15 -j 8 -v\n");
            return_code = 0;
            break;
        }
        case fmd_mode_convert: {
            fprintf(stderr, "[fmd:help] ---------- fmd:convert ----------\n");
            fprintf(stderr, "[fmd:help] fmd convert converts rgfa and FASTQ files to fmdg files and fmd query files.\n");
            fprintf(stderr, "[fmd:help] fmd convert supports two conversion modes, which are described below:\n");
            fprintf(stderr, "\trgfa2fmdg:   Converts an rGFA into an fmdg file. Throws away naming conventions and stable sequences. Everything is assumed to be on the forward strand.\n");
            fprintf(stderr, "\tfastq2query: Extracts the sequence in every read and returns one sequence per line.\n");
            fprintf(stderr, "[fmd:help] Parameters\n");
            fprintf(stderr, "\t--input   or -i: Optional parameter. Path to the input file in rGFA or FASTQ format. Default: stdin\n");
            fprintf(stderr, "\t--output  or -o: Optional parameter. Path to the input file in rGFA or FASTQ format. Default: stdout\n");
            fprintf(stderr, "\t--verbose or -v: Optional flag. Provides more information about the conversion process. Default: false\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd convert -i mygraph.rgfa -o -mygraph.fmdg -v\n");
            return_code = -1;
            break;
        }
        case fmd_mode_help: {
            fprintf(stderr, "[fmd:help] ---------- fmd:help ----------\n");
            fprintf(stderr, "[fmd:help] fmd help prints general information about how other programs under fmd work, then quits.\n");
            fprintf(stderr, "[fmd:help] To learn about a particular program, please call fmd help <program_name> through the command line.\n");
            fprintf(stderr, "[fmd:help] Any other parameters supplied are ignored.\n");
            fprintf(stderr, "[fmd:help] List of currently available programs: ");
            for(int i = 0; i < fmd_mode_no_modes; i++) {
                fprintf(stderr, "%s ", fmd_mode_names[i]);
            }
            fprintf(stderr, "\n");
            return_code = 0;
            break;
        }
        default: {
            fprintf(stderr, "[fmd:help] Program %s not recognized. Please run fmd help for more information. Quitting.\n",progname);
            return_code = -1;
            break;
        }
    }
    return return_code;
}
