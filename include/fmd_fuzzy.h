#ifndef FMD_FMD_FUZZY_H
#define FMD_FMD_FUZZY_H

#include "fmd_fmd.h"

/*
 * Warning: exponential blowup on top of exponential blowup.
 * Not for the faint of heart
 */

typedef struct fmd_fuzzy_qr_ {
    fmd_fmd_qr_t *qr;
    fmd_edit_op_t *edits;
    int_t no_edits;
} fmd_fuzzy_qr_t;

void            fmd_fuzzy_qr_init(fmd_fuzzy_qr_t **qr, fmd_string_t *query);
int             fmd_fuzzy_qr_comp(fmd_fuzzy_qr_t *q1, fmd_fuzzy_qr_t *q2);
uint_t          fmd_fuzzy_qr_hash(fmd_fuzzy_qr_t *q);
void            fmd_fuzzy_qr_free(fmd_fuzzy_qr_t *q);
fmd_fuzzy_qr_t *fmd_fuzzy_qr_copy(fmd_fuzzy_qr_t *f1);

static fmd_fstruct_t fmd_fstruct_fuzzy_qr = {
        (fcomp) fmd_fuzzy_qr_comp,
        (fhash) fmd_fuzzy_qr_hash,
        (ffree) fmd_fuzzy_qr_free,
        (fcopy) fmd_fuzzy_qr_copy,
};

typedef struct fmd_step_matcher_ {
    fmd_fmd_t *fmd;
    fmd_string_t *query;
    fmd_vector_t *fuzzy_qrs;
    fmd_vector_t *fork_nodes;
    fmd_vector_t *dead_leaves;
    int_t step;
    int_t max_edits;
} fmd_step_matcher_t;

void fmd_step_matcher_init(fmd_step_matcher_t **step, fmd_fmd_t *fmd, fmd_string_t query, int_t max_edits);
int fmd_step_matcher_comp(fmd_step_matcher_t *s1, fmd_step_matcher_t *s2);
uint_t fmd_step_matcher_hash(fmd_step_matcher_t matcher);
void fmd_step_matcher_free(fmd_step_matcher_t *step);
fmd_step_matcher_t *fmd_step_matcher_copy(fmd_step_matcher_t *step);

static fmd_fstruct_t fmd_fstruct_step_matcher = {
        (fcomp) fmd_step_matcher_comp,
        (fhash) fmd_step_matcher_hash,
        (ffree) fmd_step_matcher_free,
        (fcopy) fmd_step_matcher_copy,
};

void fmd_step_matcher_seed(fmd_step_matcher_t *step);
void fmd_step_matcher_step(fmd_step_matcher_t *step);

#endif //FMD_FMD_FUZZY_H
