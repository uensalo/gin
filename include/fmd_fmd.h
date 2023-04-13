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


typedef struct fmd_fmd_ {
    char_t c_0; // character marking the beginning of a vertex
    char_t c_1; // character marking the end of a vertex
    fmd_fmi_t *graph_fmi; // fm-index of the graph encoding
    fmd_imt_t *r2r_tree;  // translates sa ranges to sa ranges of incoming nodes
    fmd_vector_t *permutation; // permutation to make sa ranges as consecutive as possible
} fmd_fmd_t;

void fmd_fmd_init(fmd_fmd_t** fmd, fmd_graph_t *graph, fmd_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate);
void fmd_fmd_free(fmd_fmd_t *fmd);

typedef struct fmd_fmd_qr_{
    fmd_fmi_qr_t fmi_state;
    fmd_vector_t *visited;
} fmd_fmd_qr_t;

count_t fmd_fmd_query_count(fmd_fmd_t *fmd, fmd_string_t *string);
fmd_vector_t *fmd_fmd_query_locate(fmd_fmd_t *fmd, fmd_string_t *string);

fmd_vector_t *fmd_fmd_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords);

fmd_vector_t *fmd_fmd_extract_constraint_sets(fmd_graph_t *graph); // big todo :-)
                                                                   // constraint sets lead to an NP-Complete problem


#endif //FMD_FMD_FMD_H
