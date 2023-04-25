#ifndef FMD_FMD_INTERVAL_MERGE_TREE_H
#define FMD_FMD_INTERVAL_MERGE_TREE_H
#include "fmd_vector.h"

typedef struct fmd_imt_interval_{
    int_t lo;
    int_t hi;
} fmd_imt_interval_t;

void fmd_imt_interval_init(fmd_imt_interval_t **i, int_t lo, int_t hi);
void fmd_imt_interval_free(fmd_imt_interval_t *i);
fmd_imt_interval_t *fmd_imt_interval_copy(fmd_imt_interval_t *i);
uint_t fmd_imt_interval_hash(fmd_imt_interval_t *i);
int fmd_imt_interval_comp(fmd_imt_interval_t *i1, fmd_imt_interval_t *i2);

static fmd_fstruct_t fmd_fstruct_imt_interval = {
        (fcomp) fmd_imt_interval_comp,
        (fhash) fmd_imt_interval_hash,
        (ffree) fmd_imt_interval_free,
        (fcopy) fmd_imt_interval_copy
};

typedef struct fmd_imt_node_ {
    int_t lo;
    int_t hi;
    fmd_vector_t *intervals;
    struct fmd_imt_node_t *left;
    struct fmd_imt_node_t *right;
} fmd_imt_node_t;

void fmd_imt_node_init(fmd_imt_node_t **i, int_t lo, int_t hi, fmd_vector_t *intervals, fmd_imt_node_t *left, fmd_imt_node_t *right);
void fmd_imt_node_free(fmd_imt_node_t *i);
fmd_imt_node_t *fmd_imt_node_copy(fmd_imt_node_t *i);
uint_t fmd_imt_node_hash(fmd_imt_node_t *i);
int fmd_imt_node_comp(fmd_imt_node_t *i1, fmd_imt_node_t *i2);

static fmd_fstruct_t fmd_fstruct_imt_node = {
        (fcomp) fmd_imt_node_comp,
        (fhash) fmd_imt_node_hash,
        (ffree) fmd_imt_node_free,
        (fcopy) fmd_imt_node_copy
};

typedef struct fmd_imt_ {
    int_t no_keys;
    fmd_imt_node_t *root;
} fmd_imt_t;

void fmd_imt_init(fmd_imt_t **i, int_t no_keys, fmd_vector_t *kv_interval_pairs);
void fmd_imt_free(fmd_imt_t *i);
fmd_imt_t *fmd_imt_copy(fmd_imt_t *i);
uint_t fmd_imt_hash(fmd_imt_t *i);
int fmd_imt_comp(fmd_imt_t *i1, fmd_imt_t *i2);

void fmd_imt_query(fmd_imt_t *i, int_t start, int_t end, fmd_vector_t **intervals);
fmd_vector_t *fmd_imt_query_helper(fmd_imt_node_t *node, int_t lo, int_t hi);
fmd_imt_node_t *fmd_imt_init_helper(int_t lo, int_t hi, fmd_vector_t *kv_interval_pairs);
void fmd_imt_free_helper(fmd_imt_node_t *node);
fmd_vector_t *fmd_imt_merge_intervals(fmd_vector_t *i1, fmd_vector_t *i2);

fmd_imt_node_t *fmd_imt_copy_helper(fmd_imt_node_t *i);
uint_t fmd_imt_hash_helper(fmd_imt_node_t *i);
int fmd_imt_comp_helper(fmd_imt_node_t *i1, fmd_imt_node_t *i2);

static fmd_fstruct_t fmd_fstruct_imt = {
        (fcomp) fmd_imt_comp,
        (fhash) fmd_imt_hash,
        (ffree) fmd_imt_free,
        (fcopy) fmd_imt_copy
};

#endif //FMD_FMD_INTERVAL_MERGE_TREE_H
