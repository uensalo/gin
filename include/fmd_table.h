#ifndef FMD_FMD_TABLE_H
#define FMD_FMD_TABLE_H

#include "fmd_common.h"

typedef struct fmd_buffer_ {
    byte_t *buf;
    pos_t size; // number of items
    pos_t cap;  // number capacity in the number of items
} fmd_buffer_t;

typedef struct fmd_table_ {
    byte_t **buffers;
    pos_t no_buffers;
} fmd_table_t;

#endif //FMD_FMD_TABLE_H
