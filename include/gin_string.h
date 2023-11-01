/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_string.h is part of gin
 *
 * gin is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gin is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GIN_GIN_STRING_H
#define GIN_GIN_STRING_H

#include "gin_common.h"

#define GIN_STRING_INIT_SIZE 8
#define GIN_STRING_TERMINATOR '\0'
typedef char char_t;

typedef struct gin_string_ {
    char_t *seq;
    int_t size;
    int_t capacity;
} gin_string_t;

typedef enum gin_edit_op_ {
    MATCH=0,
    MISMATCH=1,
    INSERTION=2,
    DELETION=3
} gin_edit_op_t;

void gin_string_init(gin_string_t **str, int_t size);
void gin_string_init_cstr(gin_string_t **str, char *string);
void gin_string_append(gin_string_t *str, char_t c);
void gin_string_delete(gin_string_t *str, int_t pos);
void gin_string_insert(gin_string_t *str, int_t pos, char_t c);
void gin_string_substring(gin_string_t *s, int_t start, int_t end, gin_string_t **subs);
void gin_string_edit_distance(gin_string_t *str1, gin_string_t *str2, int_t *dist, gin_edit_op_t **edits, int_t *edits_len);
void gin_string_phase(gin_string_t *str1, gin_string_t *str2, gin_edit_op_t *edits, int_t edits_len, gin_string_t **p1, gin_string_t **p2);
void gin_string_phase_cigar(gin_string_t *str1, gin_string_t *str2, gin_string_t *cigar, gin_string_t **p1, gin_string_t **p2);
void gin_string_cigar(gin_edit_op_t *edits, int_t edits_len, gin_string_t **cigar);
void gin_string_free(gin_string_t *str);

gin_string_t *gin_string_copy(gin_string_t *str);
int gin_string_comp(gin_string_t *s1, gin_string_t *s2);
uint_t gin_string_hash(gin_string_t *s);

void gin_string_concat(gin_string_t **concat, gin_string_t *s1, gin_string_t *s2);
void gin_string_concat_mut(gin_string_t *s1, gin_string_t *s2);

void gin_string_kmp_lps(gin_string_t *pattern, int_t **lps);
void gin_string_kmp_search(gin_string_t *text, gin_string_t *pattern, int_t *lps, int_t **pos, int_t *num_occ);

static gin_fstruct_t gin_fstruct_string = {
        (fcomp) gin_string_comp,
        (fhash) gin_string_hash,
        (ffree) gin_string_free,
        (fcopy) gin_string_copy
};


#endif //GIN_GIN_STRING_H
