#ifndef FMD_FMD_MATCH_CHAIN_H
#define FMD_FMD_MATCH_CHAIN_H

#include "fmd_common.h"
#include "fmd_string.h"

typedef struct fmd_fmd_match_node_ {
    struct fmd_fmd_match_node_ *next;
    struct fmd_fmd_match_node_ *prev;
    fmd_string_t *matching_substring;
    int_t v_lo;
    int_t v_hi;
    int_t sa_lo;
    int_t sa_hi;
} fmd_fmd_match_node_t;
typedef struct fmd_match_chain_ {
    fmd_fmd_match_node_t *dummy; // dummy->next = head, dummy->prev = tail
    fmd_fmd_match_node_t *head;
    fmd_fmd_match_node_t *tail;
    int_t size;
} fmd_match_chain_t;
void fmd_fmd_match_chain_init(fmd_match_chain_t **list);
void fmd_fmd_match_chain_free(fmd_match_chain_t *list);
fmd_match_chain_t *fmd_fmd_match_chain_copy(fmd_match_chain_t *list);
uint_t fmd_fmd_match_chain_hash(fmd_match_chain_t *list);
int fmd_fmd_match_chain_comp(fmd_match_chain_t *l1, fmd_match_chain_t *l2);

static fmd_fstruct_t fmd_fstruct_match_list = {
        (fcomp) fmd_fmd_match_chain_comp,
        (fhash) fmd_fmd_match_chain_hash,
        (ffree) fmd_fmd_match_chain_free,
        (fcopy) fmd_fmd_match_chain_copy,
};

void fmd_fmd_match_chain_append(fmd_match_chain_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi);
void fmd_fmd_match_chain_prepend(fmd_match_chain_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi);



#endif //FMD_FMD_MATCH_CHAIN_H
