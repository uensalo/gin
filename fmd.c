#include "fmd_fmd.h"
#include <getopt.h>
#include <time.h>
#include "rgfa_parser.h"
#include "permutation_parser.h"

#define FMD_MAIN_ISA_SAMPLE_RATE_DEFAULT 256
#define FMD_MAIN_RANK_SAMPLE_RATE_DEFAULT 256
#define to_sec(t1,t2) (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1e-9
char *fmd_version = "1.0";
char* fmd_mode_names[] = {"index","query","permutation","validate","help"};

typedef enum fmd_mode_ {
    fmd_mode_index=0,
    fmd_mode_query=1,
    fmd_mode_permutation=2,
    fmd_mode_validate=3,
    fmd_mode_help=4,
    fmd_mode_no_modes=5
} fmd_mode_t;

int fmd_main_index(int argc, char **argv);
int fmd_main_query(int argc, char **argv);
int fmd_main_permutation(int argc, char **argv);
int fmd_main_validate(int argc, char **argv);
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
            return_code = fmd_main_query(argcp, argvp);
            break;
        }
        case fmd_mode_permutation : {
            return_code = fmd_main_permutation(argcp, argvp);
            break;
        }
        case fmd_mode_validate: {
            return_code = fmd_main_validate(argcp, argvp);
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
                fperm = fopen(foutput_path, "r");
                if(!fperm) {
                    fprintf(stderr, "[fmd:index] Permutation path %s could not be opened, quitting.\n", fperm_path);
                    return_code = -1;
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
            if(finput != stdin) fclose(finput);
            if(foutput != stdout) fclose(foutput);
            if(fperm) fclose(fperm);
            return_code = -1;
            return return_code;
        }
        graph = rgfa_to_fmd_graph(rgfa);
        rgfa_free(rgfa);

    } else {
        if(finput != stdin) fclose(finput);
        if(foutput != stdout) fclose(foutput);
        if(fperm) fclose(fperm);
        fprintf(stderr, "[fmd:index] fmdg parsing mode not yet implemented, quitting\n");
        return_code = -1;
        return return_code;
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
            if(foutput != stdout) fclose(foutput);
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
        int_t in_degree = 0;
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

    if(foutput != stdout) {
        foutput = fopen(foutput_path, "w");
        if(!foutput) {
            fprintf(stderr, "[fmd:index] Output path %s could not be opened, quitting.\n", foutput_path);
            free(fmd_buf);
            return_code = -1;
            return return_code;
        }
    }

    fwrite(fmd_buf,sizeof(unsigned char),fmd_buf_size,foutput);
    if(foutput != stdout) fclose(foutput);
    clock_gettime(CLOCK_REALTIME, &t2);
    write_time = to_sec(t1,t2);
    if(verbose) {
        fprintf(stderr, "[fmd:index] Resulting index size: %lld bytes (%.3lf GB).\n", fmd_buf_size, (double) fmd_buf_size * 1e-9);
        fprintf(stderr, "[fmd:index] Timings in seconds: \n");
        fprintf(stderr, "\t Parsing:  %lf", parse_time);
        fprintf(stderr, "\t Indexing: %lf", index_time);
        fprintf(stderr, "\t Writing:  %lf", write_time);
        fprintf(stderr, "[fmd:index] Total relevant runtime: %lf seconds.\n", parse_time + index_time + write_time);
    }

    free(fmd_buf);

    return return_code;
}

int fmd_main_query(int argc, char **argv) {

    return 1;
}

int fmd_main_permutation(int argc, char **argv) {

    return 1;
}

int fmd_main_validate(int argc, char **argv) {

    return 1;
}

int fmd_main_help(fmd_mode_t progmode, char *progname) {
    int return_code = 0;
    if(!progname) {
        fprintf(stderr, "%s%s%s\n%s%s",
                "[fmd:help] fmd! FM-Index like graph indexing algorithm toolkit \n",
                "[fmd:help] Needle in a haystack? More like string in a graph. Version ",
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
            fprintf(stderr, "\t--isa_sample-rate  or -s: Optional parameter. Sampling rate of the suffix array. Reducing this parameter increases query speeds at the cost of larger index files. Default = 256\n");
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
            fprintf(stderr, "[fmd:help] fmd supports many querying modes, which are described below:\n");
            fprintf(stderr, "\tcount:     Counts the occurrences of the string in the space of strings induced by the string graph. Results are returned as a single integer per line.\n");
            fprintf(stderr, "\tlocate:    Locates the vertex IDs and positions into vertices of the string matches. Results are returned as (vertex_id, index) tuples per line.\n");
            fprintf(stderr, "\tenumerate: Enumerates the paths to which vertices match. Results are returned as (vertex_id, index) (vertex_id) ... (vertex_id) (vertex_id index) per match per line.\n");
            fprintf(stderr, "[fmd:help] Parameters:\n");
            fprintf(stderr, "\t--reference or -r: Required parameter. Path to the index file. See fmd index for more help.\n");
            fprintf(stderr, "\t--input     or -i: Optional parameter. Path to the input file containing string queries, with one string per line. Default: stdin\n");
            fprintf(stderr, "\t--fastq     or -f: Optional flag.      Specifies if queries are contained in fastq format. Default: False\n");
            fprintf(stderr, "\t--output    or -o: Optional parameter. Path to the output file, produced in one of the query mode formats described above. Default: stdout\n");
            fprintf(stderr, "\t--threads   or -j: Optional parameter. Number of threads to be used for parallel querying. Default: 1\n");
            fprintf(stderr, "\t--verbose   or -v: Optional parameter. Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd query enumerate -r myindex.fmdi -i queries.fastq -f -o results.txt -j 8 -v\n");
            return_code = 0;
            break;
        }
        case fmd_mode_permutation : {
            fprintf(stderr, "[fmd:help] ---------- fmd:permutation ----------\n");
            fprintf(stderr, "[fmd:help] fmd permutation approximates a permutation of vertex indices such that they occur consecutively in the Burrows-Wheeler order.\n");
            fprintf(stderr, "[fmd:help] Computing such a permutation is NP-Complete; this program essentially runs simulated annealing and outputs a permutation.\n");
            fprintf(stderr, "[fmd:help] Parameters:\n");
            fprintf(stderr, "\t--input       or -i: Optional parameter. Path to the input file in rGFA or fmdg format. Default: stdin\n");
            fprintf(stderr, "\t--output      or -o: Optional parameter. Path to the output file, one integer as vid per line. Default: stdout\n");
            fprintf(stderr, "\t--permutation or -p: Optional parameter. Path to initial permutation file to start optimizing from. Default: Identity permutation\n");
            fprintf(stderr, "\t--time        or -t: Optional parameter. Time after which optimization is terminated in seconds. Default: 15 seconds\n");
            fprintf(stderr, "\t--update      or -u: Optional parameter. Time interval of informative prints in seconds. Default: 3 seconds\n");
            fprintf(stderr, "\t--threads     or -j: Optional parameter. Number of threads to be used for parallel cost computation. Default: 1\n");
            fprintf(stderr, "\t--verbose     or -v: Optional parameter. Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd permutation -i mygraph.fmdg -o mygraph-perm.txt -t 300 -u 15 -j 8 -v\n");
            return_code = 0;
            break;
        }
        case fmd_mode_validate: {
            fprintf(stderr, "[fmd:help] ---------- fmd:validate ----------\n");
            fprintf(stderr, "[fmd:help] Not yet implemented. Just a dud.\n");
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
