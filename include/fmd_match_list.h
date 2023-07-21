#ifndef FMD_FMD_MATCH_LIST_H
#define FMD_FMD_MATCH_LIST_H

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
typedef struct fmd_fmd_match_list_ {
    fmd_fmd_match_node_t *dummy; // dummy->next = head, dummy->prev = tail
    fmd_fmd_match_node_t *head;
    fmd_fmd_match_node_t *tail;
    int_t size;
} fmd_fmd_match_list_t;
void fmd_fmd_match_list_init(fmd_fmd_match_list_t **list);
void fmd_fmd_match_list_free(fmd_fmd_match_list_t *list);
fmd_fmd_match_list_t *fmd_fmd_match_list_copy(fmd_fmd_match_list_t *list);
uint_t fmd_fmd_match_list_hash(fmd_fmd_match_list_t *list);
int fmd_fmd_match_list_comp(fmd_fmd_match_list_t *l1, fmd_fmd_match_list_t *l2);

static fmd_fstruct_t fmd_fstruct_match_list = {
        (fcomp) fmd_fmd_match_list_comp,
        (fhash) fmd_fmd_match_list_hash,
        (ffree) fmd_fmd_match_list_free,
        (fcopy) fmd_fmd_match_list_copy,
};

void fmd_fmd_match_list_append(fmd_fmd_match_list_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi);
void fmd_fmd_match_list_prepend(fmd_fmd_match_list_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi);



#endif //FMD_FMD_MATCH_LIST_H
