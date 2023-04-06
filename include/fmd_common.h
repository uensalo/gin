#ifndef FMD_FMD_COMMON_H
#define FMD_FMD_COMMON_H
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MIN3(a,b,c) (MIN2((MIN2((a),(b))),(c)))

#define MAX2(a,b) ((a) >= (b) ? (a) : (b))
#define MAX3(a,b,c) (MAX2((MAX2((a),(b))),(c)))
#define MAX4(a,b,c,d) (MAX2((MAX3((a),(b),(c))),(d)))

typedef uint8_t  byte_t;
typedef int64_t  pos_t;
typedef uint64_t upos_t;

typedef int    (*fcomp)(void*, void*);
typedef upos_t (*fhash)(void*);
typedef void   (*ffree)(void*);
typedef void*  (*fcopy)(void*);
typedef void   (*ftrav_kv)(void *key, void *value, void *p);

typedef struct fmd_fstruct {
    fcomp comp_f;
    fhash hash_f;
    ffree free_f;
    fcopy copy_f;
} fmd_fstruct_t;

// basic functions
static int prm_comp_f(void *a, void *b) {
    return a < b ? -1 : (a == b ? 0 : 1);
}
static upos_t prm_hash_f(void *a) {
    return (pos_t)a;
}
static void prm_free(void *a) {
}
static void* prm_copy(void *val) {
    return val;
}
static fmd_fstruct_t prm_fstruct = {prm_comp_f, prm_hash_f, prm_free, prm_copy};

// swap utility
static inline void fmd_swap(void **a, void **b) {
    void *temp = *a;
    *a = *b;
    *b = temp;
}


#endif //FMD_FMD_COMMON_H
