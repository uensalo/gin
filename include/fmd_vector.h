#ifndef FMD_FMD_VECTOR_H
#define FMD_FMD_VECTOR_H
#include "fmd_common.h"

#define FMD_VECTOR_INIT_SIZE 8

typedef struct fmd_vector_ {
    void **data;
    pos_t size;
    pos_t capacity;
} fmd_vector_t;

#endif //FMD_FMD_VECTOR_H
