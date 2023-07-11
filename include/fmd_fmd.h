#ifndef FMD_FMD_FMD_H
#define FMD_FMD_FMD_H
#include "fmd_common.h"
#include "fmd_graph.h"
#include "fmd_interval_merge_tree.h"
#include "fmd_fmi.h"
#include "fmd_table.h"
#include "fmd_vector.h"
#include "assert.h"

// separator characters
#define FMD_FMD_DEFAULT_c_0 '('
#define FMD_FMD_DEFAULT_c_1 ')'
// permutation encoding characters
#define FMD_FMD_DEFAULT_a_0 ','
#define FMD_FMD_DEFAULT_a_1 '.'
// bit field lengths for serialization
#define FMD_FMD_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_SPECIAL_CHAR_BIT_LENGTH 40
#define FMD_FMD_NO_VERTICES_BIT_LENGTH 40
#define FMD_FMD_PERMUTATION_BIT_LENGTH 40
#define FMD_FMD_BWT_TO_VID_BIT_LENGTH 40
#define FMD_FMD_FMI_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_IMT_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH 32
#define FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH 40

typedef struct fmd_fork_node_{
    struct fmd_fork_node_t *parent;
    int_t vertex_lo;
    int_t vertex_hi;
    int_t sa_lo;
    int_t sa_hi;
    int_t pos;
    bool is_leaf;
    bool is_dead;
    bool is_cadet; // used to be called is bastard
} fmd_fork_node_t;
fmd_fork_node_t *fmd_fork_node_init(fmd_fork_node_t *parent, int_t vlo, int_t vhi, int_t pos, bool is_leaf, bool is_dead, bool is_cadet);
void fmd_fork_node_free(fmd_fork_node_t *node);
fmd_fork_node_t *fmd_fork_node_copy(fmd_fork_node_t *node);
uint_t fmd_fork_node_hash(fmd_fork_node_t *node);
int fmd_fork_node_comp(fmd_fork_node_t *n1, fmd_fork_node_t *n2);

static fmd_fstruct_t fmd_fstruct_fork_node = {
        (fcomp) fmd_fork_node_comp,
        (fhash) fmd_fork_node_hash,
        (ffree) fmd_fork_node_free,
        (fcopy) fmd_fork_node_copy,
};

typedef struct fmd_fmd_qr_{
    fmd_fork_node_t *cur_fork;
    int_t lo;
    int_t hi;
    int_t pos;
    fmd_string_t *pattern;
} fmd_fmd_qr_t;
fmd_fmd_qr_t *fmd_fmd_qr_init(fmd_fork_node_t *cur_fork, int_t lo, int_t hi, int_t pos, fmd_string_t *pattern);
void fmd_fmd_qr_free(fmd_fmd_qr_t *q);

typedef struct fmd_fmd_ {
    char_t c_0; // character marking the beginning of a vertex
    char_t c_1; // character marking the end of a vertex
    fmd_vector_t *permutation; // permutation to make sa ranges as consecutive as possible
    fmd_vector_t *bwt_to_vid; // converts c0 ranks to text ranks, i.e. vids
    fmd_fmi_t *graph_fmi; // fm-index of the graph encoding
    fmd_imt_t *r2r_tree;  // translates sa ranges to sa ranges of incoming nodes
} fmd_fmd_t;

void fmd_fmd_init(fmd_fmd_t** fmd, fmd_graph_t *graph, fmd_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate);
void fmd_fmd_free(fmd_fmd_t *fmd);

int fmd_fmd_comp(fmd_fmd_t *f1, fmd_fmd_t *f2);

count_t fmd_fmd_query_count(fmd_fmd_t *fmd, fmd_string_t *string);
fmd_vector_t *fmd_fmd_query_locate_basic(fmd_fmd_t *fmd, fmd_string_t *string);
void fmd_fmd_query_locate_paths(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_matches, fmd_vector_t **paths, fmd_vector_t **dead_ends);
void fmd_fmd_locate_paths_result_free(fmd_vector_t *paths, fmd_vector_t *dead_ends);
void fmd_fmd_query_locate_paths_stats(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **paths, fmd_vector_t **dead_ends, int_t *no_forks);

#ifdef FMD_OMP
void fmd_fmd_query_locate_paths_omp(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_matches, fmd_vector_t **paths, fmd_vector_t **dead_ends, int_t num_threads);
void fmd_fmd_query_locate_paths_process_query_record(fmd_fmd_t *fmd, fmd_fmd_qr_t *rec, int_t max_matches, fmd_vector_t *exact_matches, fmd_vector_t *partial_matches);
#endif

bool fmd_fmd_advance_query(fmd_fmi_t *fmi, fmd_fmd_qr_t *qr);
bool fmd_fmd_query_precedence_range(fmd_fmi_t *fmi, fmd_fmd_qr_t *qr, char_t c, int_t *lo, int_t *hi);

fmd_vector_t *fmd_fmd_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords);

void fmd_fmd_serialize_from_buffer(fmd_fmd_t **fmd_ret, unsigned char *buf, uint64_t buf_size);
/**
 * Write order to buffer:
 *  - FMD_FMD_FMI_NO_BITS_BIT_LENGTH : fmd_bit_length // Number of bits for the entire data structure
 *  - FMD_FMD_SPECIAL_CHAR_BIT_LENGTH : c_0 // Special character char c_0
 *  - FMD_FMD_SPECIAL_CHAR_BIT_LENGTH : c_1 // Special character char c_1
 *  - FMD_FMD_NO_VERTICES_BIT_LENGTH : no_vertices // number of vertices in the graph of the fmd
 *  - for i = 0 to the number of vertices:
 *      - FMD_FMD_PERMUTATION_BIT_LENGTH permutation[i] // Permutation entry
 *  - for i = 0 to the number of vertices:
 *      - FMD_FMD_BWT_TO_VID_BIT_LENGTH bwt_to_vid[i] // Rank translation entry
 *  - FMD_FMD_FMI_NO_BITS_BIT_LENGTH : fmi_bit_length // bit-length of the FMI bitstream
 *  - fmi_bit_length : fmi // the buffer for the fm-index, word aligned
 *  - FMD_FMD_IMT_NO_BITS_BIT_LENGTH : imt_bit_length // bit_length of the IMT buffer
 *  - imt_bit_length: imt // imt buffer
 * @param fmd
 * @param buf
 * @param buf_size
 */
void fmd_fmd_serialize_to_buffer(fmd_fmd_t *fmd, unsigned char **buf_ret, uint64_t *buf_size_re);
void fmd_fmd_serialize_to_buffer_imt_helper(fmd_imt_node_t *node, fmd_bs_t *bs, uint_t *widx);


#endif //FMD_FMD_FMD_H
