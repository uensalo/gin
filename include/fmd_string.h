#ifndef FMD_FMD_STRING_H
#define FMD_FMD_STRING_H

#include "fmd_common.h"

#define FMD_STRING_INIT_SIZE 8
#define FMD_STRING_TERMINATOR '\0'
typedef char char_t;

typedef struct fmd_string_ {
    char_t *seq;
    int_t size;
    int_t capacity;
} fmd_string_t;

typedef enum fmd_edit_op_ {
    MATCH=0,
    MISMATCH=1,
    INSERTION=2,
    DELETION=3
} fmd_edit_op_t;

void fmd_string_init(fmd_string_t **str, int_t size);
void fmd_string_init_cstr(fmd_string_t **str, char *string);
void fmd_string_append(fmd_string_t *str, char_t c);
void fmd_string_delete(fmd_string_t *str, int_t pos);
void fmd_string_insert(fmd_string_t *str, int_t pos, char_t c);
void fmd_string_substring(fmd_string_t *s, int_t start, int_t end, fmd_string_t **subs);
void fmd_string_edit_distance(fmd_string_t *str1, fmd_string_t *str2, int_t *dist, fmd_edit_op_t **edits, int_t *edits_len);
void fmd_string_phase(fmd_string_t *str1, fmd_string_t *str2, fmd_edit_op_t *edits, int_t edits_len, fmd_string_t **p1, fmd_string_t **p2);
void fmd_string_phase_cigar(fmd_string_t *str1, fmd_string_t *str2, fmd_string_t *cigar, fmd_string_t **p1, fmd_string_t **p2);
void fmd_string_cigar(fmd_edit_op_t *edits, int_t edits_len, fmd_string_t **cigar);
void fmd_string_free(fmd_string_t *str);

fmd_string_t *fmd_string_copy(fmd_string_t *str);
int fmd_string_comp(fmd_string_t *s1, fmd_string_t *s2);
uint_t fmd_string_hash(fmd_string_t *s);

void fmd_string_concat(fmd_string_t **concat, fmd_string_t *s1, fmd_string_t *s2);
void fmd_string_concat_mut(fmd_string_t *s1, fmd_string_t *s2);

void fmd_string_kmp_lps(fmd_string_t *pattern, int_t **lps);
void fmd_string_kmp_search(fmd_string_t *text, fmd_string_t *pattern, int_t *lps, int_t **pos, int_t *num_occ);

static fmd_fstruct_t fmd_fstruct_string = {
        (fcomp) fmd_string_comp,
        (fhash) fmd_string_hash,
        (ffree) fmd_string_free,
        (fcopy) fmd_string_copy
};


#endif //FMD_FMD_STRING_H
