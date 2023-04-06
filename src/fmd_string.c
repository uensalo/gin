#include "fmd_string.h"

void fmd_string_init(fmd_string_t **str, pos_t size) {
    fmd_string_t *s = calloc(1, sizeof(fmd_string_t));
    if(!s) return;
    s->seq = calloc((s->capacity=size)+1, sizeof(fmd_char_t));
    *str = s;
}

void fmd_string_init_cstr(fmd_string_t **str, char *string) {
    size_t len = strlen(string);
    fmd_string_t *s = calloc(1, sizeof(fmd_string_t));
    s->seq = calloc((s->capacity=2*(pos_t)len)+1, sizeof(fmd_char_t));
    s->size = (pos_t)len;
    memcpy(s->seq, string, s->size * sizeof(fmd_char_t));
    *str = s;
}

void fmd_string_append(fmd_string_t *str, fmd_char_t c) {
    if(str->size == str->capacity) {
        str->seq = realloc(str->seq, ((str->capacity *= 2)+1) * sizeof(fmd_char_t));
        memset(str->seq + str->size, 0, (str->size+1) * sizeof(fmd_char_t));
    }
    str->seq[str->size++] = c;
}

void fmd_string_delete(fmd_string_t *str, pos_t pos) {
    for(pos_t i = pos; i < str->size-1; i++) {
        str->seq[i] = str->seq[i+1];
    }
    str->size--;
    str->seq[str->size] = 0;
}

void fmd_string_insert(fmd_string_t *str, pos_t pos, fmd_char_t c) {
    if(str->size == str->capacity) {
        str->seq = realloc(str->seq, ((str->capacity *= 2)+1) * sizeof(fmd_char_t));
        memset(str->seq + str->size, 0, (str->size+1) * sizeof(fmd_char_t));
    }
    for(pos_t i = str->size; i > pos; i--) {
        str->seq[i] = str->seq[i-1];
    }
    str->seq[pos] = c;
    str->size++;
    str->seq[str->size] = 0;
}

void fmd_string_substring(fmd_string_t *s, pos_t start, pos_t end, fmd_string_t **subs) {
    if (start >= end || start >= s->size || end > s->size) {
        *subs = NULL;
        return;
    }
    pos_t substring_length = end - start;
    fmd_string_init(subs, substring_length);
    memcpy((*subs)->seq, s->seq + start, substring_length * sizeof(fmd_char_t));
    (*subs)->size = substring_length;
    (*subs)->seq[substring_length] = 0;
}

void fmd_string_edit_distance(fmd_string_t *str1, fmd_string_t *str2, pos_t *dist, fmd_edit_op_t **edits, pos_t *edits_len) {
    // init matrix
#define VAL(m, i ,j, n) (m)[(i) * (n) + (j)]
    pos_t n = str1->size + 1;
    pos_t m = str2->size + 1;
    pos_t *dp = calloc(n*m, sizeof(pos_t));

    fmd_edit_op_t *cigar = calloc(n+m, sizeof(fmd_edit_op_t));

    for(pos_t i = 1; i < n; i++) {
        VAL(dp, i, 0, m) = i;
    }
    for(pos_t j = 1; j < m; j++) {
        VAL(dp, 0, j, m) = j;
    }
    for(pos_t i = 1; i < n; i++) {
        for(pos_t j = 1; j < m; j++) {
            if(str1->seq[i-1] == str2->seq[j-1]) {
                VAL(dp, i, j, m) = VAL(dp, i-1,j-1, m);
            } else {
                pos_t mismatch = VAL(dp, i-1, j-1, m) + 1;
                pos_t insert = VAL(dp, i-1, j, m) + 1;
                pos_t delete = VAL(dp, i, j-1, m) + 1;
                pos_t minval =  MIN3(mismatch, insert, delete);
                VAL(dp, i, j, m) = minval;
            }
        }
    }
    // retval and traceback
    pos_t i = n-1;
    pos_t j = m-1;
    pos_t t = 0;
    while(i > 0 && j > 0) {
        pos_t overwrite = VAL(dp, i-1, j-1, m);
        pos_t insert = VAL(dp, i-1, j, m);
        pos_t delete = VAL(dp, i, j-1, m);
        pos_t minval = MIN3(overwrite, insert, delete);
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
    pos_t a = t-1;
    pos_t b = 0;
    while(a > b) {
        fmd_edit_op_t tmp = cigar[a];
        cigar[a] = cigar[b];
        cigar[b] = tmp;
        --a;++b;
    }
    free(dp);
#undef VAL
}

void fmd_string_phase(fmd_string_t *str1, fmd_string_t *str2, fmd_edit_op_t *edits, pos_t edits_len, fmd_string_t **p1, fmd_string_t **p2) {
    fmd_string_t *s1,*s2;
    fmd_string_init(&s1, edits_len+1);
    fmd_string_init(&s2, edits_len+1);

    pos_t r1 = 0;
    pos_t r2 = 0;
    for(pos_t i = 0; i < edits_len; i++) {
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

    pos_t r1 = 0;
    pos_t r2 = 0;
    for(pos_t i = 0; i < cigar->size; i++) {
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

void fmd_string_cigar(fmd_edit_op_t *edits, pos_t edits_len, fmd_string_t **cigar) {
    fmd_string_t *out;
    fmd_string_init(&out, edits_len);
    for (pos_t i = 0; i < edits_len; i++) {
        fmd_char_t c;
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
    copy->seq = calloc(str->capacity, sizeof(fmd_char_t));
    if(!copy->seq) {
        free(copy);
        return NULL;
    }
    memcpy(copy->seq, str->seq, sizeof(fmd_char_t) * str->size);
    copy->size = str->size;
    copy->capacity = str->capacity;
    return copy;
}

int fmd_string_comp(fmd_string_t *s1, fmd_string_t *s2) {
    pos_t min_len = s1->size > s2->size ? s2->size : s1->size;
    pos_t i;
    for(i = 0; i < min_len; i++) {
        if(s1->seq[i] < s2->seq[i])
            return -1;
        else if(s1->seq[i] > s2->seq[i])
            return 1;
    }
    if (s1->size == s2->size) return 0;
    return s1->size > s2->size ? -1 : 1;
}

upos_t fmd_string_hash(fmd_string_t *s) {
    const upos_t prime = 1099511628211LLU;
    upos_t hash = 14695981039346656037LLU;
    uint64_t *seq64 = (uint64_t*)s->seq;
    int block_size = sizeof(uint64_t) / sizeof(fmd_char_t);
    for (pos_t i = 0; i < s->size/block_size; ++i) {
        hash ^= seq64[i];
        hash *= prime;
    }
    for(pos_t i = (s->size/block_size)*block_size; s->size; i++) {
        hash ^= s->seq[i];
        hash *= prime;
    }
    return hash;
}