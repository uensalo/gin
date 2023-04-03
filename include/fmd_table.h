#ifndef FMD_FMD_TABLE_H
#define FMD_FMD_TABLE_H

#include "fmd_common.h"
#include "fmd_tree.h"

#define FMD_HT_INIT_SIZE 53
#define FMD_REHASH_FACTOR 2

typedef struct fmd_table_ {
    pos_t size;
    pos_t capacity;
    fmd_fstruct_t *key_f;
    fmd_fstruct_t *val_f;
    fmd_tree_node_t **roots;
    pos_t *items_per_bucket;
} fmd_table_t;

unsigned long djb2(const char *str) {
    unsigned long hash = 5381;
    int c;

    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c;
    }

    return hash;
}


#endif //FMD_FMD_TABLE_H
