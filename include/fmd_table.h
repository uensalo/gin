#ifndef FMD_FMD_TABLE_H
#define FMD_FMD_TABLE_H

#include "fmd_common.h"

typedef int (*key_cmp)(void*, void*);
typedef pos_t (*key_hash)(void*);

typedef struct fmd_kv_ {
     void *k;
     void *v;
} fmd_kv_t;

typedef struct fmd_buffer_ {
    fmd_kv_t *buf;
    pos_t size; // number of items
    pos_t cap;  // number capacity in the number of items
} fmd_vec_t;

typedef struct fmd_table_ {
    fmd_vec_t *buckets;
    pos_t no_buckets;
    key_cmp cmp;
    key_hash hsh;
    pos_t no_items;
} fmd_table_t;


#endif //FMD_FMD_TABLE_H
