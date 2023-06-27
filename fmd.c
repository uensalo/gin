#include "fmd_fmd.h"
#include <getopt.h>

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


    return 1;
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
            fprintf(stderr, "\t--output           or -o: Optional parameter. Path to the output file, produced in binary fmdi format. Default: stdout\n");
            fprintf(stderr, "\t--permutation      or -p: Optional parameter. Path to the permutation file. See fmd permutation for more help. Default = Identity permutation\n");
            fprintf(stderr, "\t--isa_sample-rate  or -s: Optional parameter. Sampling rate of the suffix array. Reducing this parameter increases query speeds at the cost of larger index files. Default = 256\n");
            fprintf(stderr, "\t--rank-sample-rate or -r: Optional parameter. Frequency of rank caches. Reducing this parameter increases query speeds at the cost of larger index files. Default = 256\n");
            fprintf(stderr, "\t--verbose          or -v: Optional flag.      Provides more information (time, progress, memory requirements) about the indexing process.\n");
            fprintf(stderr, "[fmd:help] Example invocation: fmd index -i mygraph.fmdg -o mygraph.fmdi -p myperm -s 64 -r 64 -v\n");
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
