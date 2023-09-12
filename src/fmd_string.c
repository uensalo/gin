/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_string.c is part of fmd
 *
 * fmd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * fmd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "fmd_string.h"

void fmd_string_init(fmd_string_t **str, int_t size) {
    fmd_string_t *s = calloc(1, sizeof(fmd_string_t));
    if(!s) return;
    s->seq = calloc((s->capacity=size)+1, sizeof(char_t));
    *str = s;
}

void fmd_string_init_cstr(fmd_string_t **str, char *string) {
    size_t len = strlen(string);
    fmd_string_t *s = calloc(1, sizeof(fmd_string_t));
    s->seq = calloc((s->capacity=2*(int_t)len)+1, sizeof(char_t));
    s->size = (int_t)len;
    memcpy(s->seq, string, s->size * sizeof(char_t));
    *str = s;
}

void fmd_string_append(fmd_string_t *str, char_t c) {
    if(str->size == str->capacity) {
        str->seq = realloc(str->seq, ((str->capacity *= 2)+1) * sizeof(char_t));
        memset(str->seq + str->size, FMD_STRING_TERMINATOR, (str->size+1) * sizeof(char_t));
    }
    str->seq[str->size++] = c;
}

void fmd_string_delete(fmd_string_t *str, int_t pos) {
    for(int_t i = pos; i < str->size-1; i++) {
        str->seq[i] = str->seq[i+1];
    }
    str->size--;
    str->seq[str->size] = FMD_STRING_TERMINATOR;
}

void fmd_string_insert(fmd_string_t *str, int_t pos, char_t c) {
    if(str->size == str->capacity) {
        str->seq = realloc(str->seq, ((str->capacity *= 2)+1) * sizeof(char_t));
        memset(str->seq + str->size, 0, (str->size+1) * sizeof(char_t));
    }
    for(int_t i = str->size; i > pos; i--) {
        str->seq[i] = str->seq[i-1];
    }
    str->seq[pos] = c;
    str->size++;
    str->seq[str->size] = FMD_STRING_TERMINATOR;
}

void fmd_string_substring(fmd_string_t *s, int_t start, int_t end, fmd_string_t **subs) {
    if (start >= end || start >= s->size || end > s->size) {
        *subs = NULL;
        return;
    }
    int_t substring_length = end - start;
    fmd_string_init(subs, substring_length);
    memcpy((*subs)->seq, s->seq + start, substring_length * sizeof(char_t));
    (*subs)->size = substring_length;
    (*subs)->seq[substring_length] = 0;
}

void fmd_string_edit_distance(fmd_string_t *str1, fmd_string_t *str2, int_t *dist, fmd_edit_op_t **edits, int_t *edits_len) {
    // init matrix
#define VAL(m, i ,j, n) (m)[(i) * (n) + (j)]
    int_t n = str1->size + 1;
    int_t m = str2->size + 1;
    int_t *dp = calloc(n*m, sizeof(int_t));

    fmd_edit_op_t *cigar = calloc(n+m, sizeof(fmd_edit_op_t));

    for(int_t i = 1; i < n; i++) {
        VAL(dp, i, 0, m) = i;
    }
    for(int_t j = 1; j < m; j++) {
        VAL(dp, 0, j, m) = j;
    }
    for(int_t i = 1; i < n; i++) {
        for(int_t j = 1; j < m; j++) {
            if(str1->seq[i-1] == str2->seq[j-1]) {
                VAL(dp, i, j, m) = VAL(dp, i-1,j-1, m);
            } else {
                int_t mismatch = VAL(dp, i-1, j-1, m) + 1;
                int_t insert = VAL(dp, i-1, j, m) + 1;
                int_t delete = VAL(dp, i, j-1, m) + 1;
                int_t minval =  MIN3(mismatch, insert, delete);
                VAL(dp, i, j, m) = minval;
            }
        }
    }
    // retval and traceback
    int_t i = n-1;
    int_t j = m-1;
    int_t t = 0;
    while(i > 0 && j > 0) {
        int_t overwrite = VAL(dp, i-1, j-1, m);
        int_t insert = VAL(dp, i-1, j, m);
        int_t delete = VAL(dp, i, j-1, m);
        int_t minval = MIN3(overwrite, insert, delete);
        fmd_edit_op_t op;
        if(minval == overwrite) {
            op = str1->seq[i-1] == str2->seq[j-1] ? MATCH : MISMATCH;
            --i;--j;
        } else if(minval == insert) {
            op = INSERTION;
            --i;
        } else if(minval == delete) {
            op = DELETION;
            --j;
        }
        cigar[t++] = op;
    }
    fmd_edit_op_t op;
    while(i>0) {
        op = INSERTION;
        cigar[t++] = op;
        --i;
    }
    while(j>0) {
        op = DELETION;
        cigar[t++] = op;
        --j;
    }
    *dist = VAL(dp, n-1, m-1, m);
    *edits_len = t;
    *edits = cigar;
    // reverse the cigar string
    int_t a = t-1;
    int_t b = 0;
    while(a > b) {
        fmd_edit_op_t tmp = cigar[a];
        cigar[a] = cigar[b];
        cigar[b] = tmp;
        --a;++b;
    }
    free(dp);
#undef VAL
}

void fmd_string_phase(fmd_string_t *str1, fmd_string_t *str2, fmd_edit_op_t *edits, int_t edits_len, fmd_string_t **p1, fmd_string_t **p2) {
    fmd_string_t *s1,*s2;
    fmd_string_init(&s1, edits_len+1);
    fmd_string_init(&s2, edits_len+1);

    int_t r1 = 0;
    int_t r2 = 0;
    for(int_t i = 0; i < edits_len; i++) {
        switch (edits[i]) {
            case MATCH:
            case MISMATCH: {
                s1->seq[s1->size++]=str1->seq[r1++];
                s2->seq[s2->size++]=str2->seq[r2++];
                break;
            }
            case DELETION: {
                s1->seq[s1->size++]='-';
                s2->seq[s2->size++]=str2->seq[r2++];
                break;
            }
            case INSERTION: {
                s1->seq[s1->size++]=str1->seq[r1++];
                s2->seq[s2->size++]='-';
                break;
            }
        }
    }
    *p1 = s1;
    *p2 = s2;
}

void fmd_string_phase_cigar(fmd_string_t *str1, fmd_string_t *str2, fmd_string_t *cigar, fmd_string_t **p1, fmd_string_t **p2) {
    fmd_string_t *s1,*s2;
    fmd_string_init(&s1, cigar->size+1);
    fmd_string_init(&s2, cigar->size+1);

    int_t r1 = 0;
    int_t r2 = 0;
    for(int_t i = 0; i < cigar->size; i++) {
        switch (cigar->seq[i]) {
            case 'M':
            case 'X': {
                s1->seq[s1->size++]=str1->seq[r1++];
                s2->seq[s2->size++]=str2->seq[r2++];
                break;
            }
            case 'D': {
                s1->seq[s1->size++]='-';
                s2->seq[s2->size++]=str2->seq[r2++];
                break;
            }
            case 'I': {
                s1->seq[s1->size++]=str1->seq[r1++];
                s2->seq[s2->size++]='-';
                break;
            }
        }
    }
    *p1 = s1;
    *p2 = s2;
}

void fmd_string_cigar(fmd_edit_op_t *edits, int_t edits_len, fmd_string_t **cigar) {
    fmd_string_t *out;
    fmd_string_init(&out, edits_len);
    for (int_t i = 0; i < edits_len; i++) {
        char_t c;
        switch (edits[i]) {
            case MATCH: {
                c = 'M';
                break;
            }
            case MISMATCH: {
                c = 'X';
                break;
            }
            case INSERTION: {
                c = 'I';
                break;
            }
            case DELETION: {
                c = 'D';
                break;
            }
        }
        fmd_string_append(out, c);
    }
    *cigar = out;
}

void fmd_string_free(fmd_string_t *str) {
    if(str) {
        free(str->seq);
        free(str);
    }
}

fmd_string_t *fmd_string_copy(fmd_string_t *str) {
    fmd_string_t *copy = calloc(1, sizeof(fmd_string_t));
    if(!copy) return NULL;
    copy->seq = calloc(str->capacity+1, sizeof(char_t));
    if(!copy->seq) {
        free(copy);
        return NULL;
    }
    memcpy(copy->seq, str->seq, sizeof(char_t) * str->size);
    copy->size = str->size;
    copy->capacity = str->capacity;
    return copy;
}

int fmd_string_comp(fmd_string_t *s1, fmd_string_t *s2) {
    if(!s1 && !s2) return 0;
    int_t min_len = s1->size > s2->size ? s2->size : s1->size;
    int_t i;
    for(i = 0; i < min_len; i++) {
        if(s1->seq[i] < s2->seq[i])
            return -1;
        else if(s1->seq[i] > s2->seq[i])
            return 1;
    }
    if (s1->size == s2->size) return 0;
    return s1->size > s2->size ? -1 : 1;
}

uint_t fmd_string_hash(fmd_string_t *s) {
    const uint_t prime = 1099511628211LLU;
    uint_t hash = 14695981039346656037LLU;
    uint64_t *seq64 = (uint64_t*)s->seq;
    int block_size = sizeof(uint64_t) / sizeof(char_t);
    for (int_t i = 0; i < s->size/block_size; ++i) {
        hash ^= seq64[i];
        hash *= prime;
    }
    for(int_t i = (s->size/block_size)*block_size; i < s->size; i++) {
        hash ^= s->seq[i];
        hash *= prime;
    }
    return hash;
}

void fmd_string_concat(fmd_string_t **concat, fmd_string_t *s1, fmd_string_t *s2) {
    int_t total_size = s1->size + s2->size;
    fmd_string_init(concat, total_size);
    memcpy((*concat)->seq, s1->seq, s1->size * sizeof(char_t));
    memcpy((*concat)->seq + s1->size, s2->seq, s2->size * sizeof(char_t));
    (*concat)->size = total_size;
    (*concat)->seq[total_size] = FMD_STRING_TERMINATOR;
}

void fmd_string_concat_mut(fmd_string_t *s1, fmd_string_t *s2) {
    int_t total_size = s1->size + s2->size;
    if (total_size > s1->capacity) {
        s1->seq = realloc(s1->seq, (total_size + 1) * sizeof(char_t));
        s1->capacity = total_size;
    }
    memcpy(s1->seq + s1->size, s2->seq, s2->size * sizeof(char_t));
    s1->size = total_size;
    s1->seq[total_size] = FMD_STRING_TERMINATOR;
}

void fmd_string_kmp_lps(fmd_string_t *pattern, int_t **lps_rval) {
    int_t *lps = calloc(pattern->size, sizeof(int_t));
    int_t length = 0; // length of the previous longest prefix suffix
    lps[0] = 0; // lps[0] is always 0

    for (int_t i = 1; i < pattern->size; i++) {
        while (length > 0 && pattern->seq[i] != pattern->seq[length]) {
            length = lps[length-1];
        }
        if (pattern->seq[i] == pattern->seq[length]) {
            length++;
        }
        lps[i] = length;
    }
    *lps_rval = lps;
}

void fmd_string_kmp_search(fmd_string_t *text, fmd_string_t *pattern, int_t *lps_in, int_t **pos_rval, int_t *occ_rval) {
    int_t *lps;
    if(!lps_in) {
        fmd_string_kmp_lps(pattern, &lps);
    } else {
        lps = lps_in;
    }
    int_t occ = 0;
    int_t cap = 8;
    int_t *pos = malloc(cap*sizeof(int_t));

    if (pattern->size == 0 || text->size == 0) {
        *pos_rval = NULL;
        *occ_rval = 0;
    }

    int_t i = 0; // text index
    int_t j = 0; // pattern index

    while (i < text->size) {
        if (pattern->seq[j] == text->seq[i]) {
            j++;
            i++;
        }
        if (j == pattern->size) {
            if(occ == cap) {
                pos = realloc(pos, (cap *= 2) * sizeof(int_t));
            }
            pos[occ++] = i - j;
            j = lps[j - 1]; // reset j to find rest of the matches
        } else if (i < text->size && pattern->seq[j] != text->seq[i]) {
            if (j != 0) {
                j = lps[j-1];
            } else {
                i++;
            }
        }
    }
    *pos_rval = realloc(pos,occ*sizeof(int_t));
    *occ_rval = occ;
    if(!lps_in) free(lps);
}

