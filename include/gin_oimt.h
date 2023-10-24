#ifndef GIN_GIN_OIMT_H
#define GIN_GIN_OIMT_H
#include "gin_interval_merge_tree.h"
#include "gin_graph.h"

typedef struct gin_oimt_node_ {
    int_t lo;
    int_t hi;
    gin_vector_t **interval_buckets;
    struct gin_oimt_node_ *left;
    struct gin_oimt_node_ *right;
} gin_oimt_node_t;

typedef struct gin_oimt_ {
    int_t no_keys; // statistics purposes
    gin_oimt_node_t *root;
    int_t *alphabet;
    int_t alphabet_size;
    int_t *c2e;
} gin_oimt_t;

void gin_oimt_init(gin_imt_t *imt_node, int_t *vertex_last_char_enc , int_t *alphabet, int_t alphabet_size, gin_oimt_t **oimt);
void gin_oimt_free(gin_oimt_t *oimt);

//gin_oimt_t* gin_oimt_copy(gin_oimt_t *o);
//uint_t gin_oimt_hash(gin_oimt_t *o);
//void gin_oimt_comp(gin_oimt_t *o1, gin_oimt_t *o2);

gin_oimt_node_t *gin_oimt_init_helper(gin_imt_node_t *imt, int_t *vertex_last_char_enc, int_t *alphabet, int_t alphabet_size);
void gin_oimt_query(gin_oimt_t *oimt, int_t start, int_t end, char_t c, int_t no_max_intervals, gin_vector_t **intervals);
void gin_oimt_query_helper(gin_oimt_node_t *node, int_t lo, int_t hi, int_t enc, int_t no_max_intervals, int_t *no_cur_intervals, gin_vector_t *merge_list);
void gin_oimt_free_helper(gin_oimt_node_t *node, int_t alphabet_size);


#endif //GIN_GIN_OIMT_H
