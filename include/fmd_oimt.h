#ifndef FMD_FMD_OIMT_H
#define FMD_FMD_OIMT_H
#include "fmd_interval_merge_tree.h"
#include "fmd_graph.h"

typedef struct fmd_oimt_node_ {
    int_t lo;
    int_t hi;
    fmd_vector_t **interval_buckets;
    struct fmd_oimt_node_ *left;
    struct fmd_oimt_node_ *right;
} fmd_oimt_node_t;

typedef struct fmd_oimt_ {
    int_t no_keys; // statistics purposes
    fmd_oimt_node_t *root;
    int_t *alphabet;
    int_t alphabet_size;
    int_t *c2e;
} fmd_oimt_t;

void fmd_oimt_init(fmd_imt_t *imt_node, int_t *vertex_last_char_enc , int_t *alphabet, int_t alphabet_size, fmd_oimt_t **oimt);
void fmd_oimt_free(fmd_oimt_t *oimt);

//fmd_oimt_t* fmd_oimt_copy(fmd_oimt_t *o);
//uint_t fmd_oimt_hash(fmd_oimt_t *o);
//void fmd_oimt_comp(fmd_oimt_t *o1, fmd_oimt_t *o2);

fmd_oimt_node_t *fmd_oimt_init_helper(fmd_imt_node_t *imt, int_t *vertex_last_char_enc, int_t *alphabet, int_t alphabet_size);
void fmd_oimt_query(fmd_oimt_t *oimt, int_t start, int_t end, char_t c, int_t no_max_intervals, fmd_vector_t **intervals);
void fmd_oimt_query_helper(fmd_oimt_node_t *node, int_t lo, int_t hi, int_t enc, int_t no_max_intervals, int_t *no_cur_intervals, fmd_vector_t *merge_list);
void fmd_oimt_free_helper(fmd_oimt_node_t *node, int_t alphabet_size);


#endif //FMD_FMD_OIMT_H
