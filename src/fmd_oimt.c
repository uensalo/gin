#include "fmd_oimt.h"

void fmd_oimt_init(fmd_imt_t *imt, int_t *vertex_last_char_enc, int_t *alphabet, int_t alphabet_size, fmd_oimt_t **oimt) {
    fmd_oimt_t *o = calloc(1, sizeof(fmd_oimt_t));
    o->root = fmd_oimt_init_helper(imt->root, vertex_last_char_enc, alphabet, alphabet_size);
    o->no_keys = imt->no_keys;
    o->alphabet = calloc(alphabet_size, sizeof(int_t));
    memcpy(o->alphabet, alphabet, alphabet_size * sizeof(int_t));
    o->alphabet_size = alphabet_size;

    int_t t = 0;
    o->c2e = calloc(FMD_MAX_ALPHABET_SIZE, sizeof(int_t));
    for(int_t i = 0; i < alphabet_size; i++) {
        o->c2e[alphabet[i]] = t;
        t++;
    }
    *oimt = o;
}

fmd_oimt_node_t *fmd_oimt_init_helper(fmd_imt_node_t *imt_node,
                                      int_t *vertex_last_char_enc,
                                      int_t *alphabet,
                                      int_t alphabet_size) {
    if(!imt_node) return NULL;
    fmd_oimt_node_t *o = calloc(1, sizeof(fmd_oimt_node_t));
    o->lo = imt_node->lo;
    o->hi = imt_node->hi;
    // partition the alphabet
    o->interval_buckets = calloc(alphabet_size, sizeof(fmd_vector_t*));
    fmd_vector_t **tmp_lists = calloc(alphabet_size, sizeof(fmd_vector_t*));
    for(int_t i = 0; i < alphabet_size; i++) {
        fmd_vector_init(&tmp_lists[i], FMD_VECTOR_INIT_SIZE, &fmd_fstruct_imt_interval);
    }

    for(int_t i = 0; i < imt_node->intervals->size; i++) {
        fmd_imt_interval_t *interval = imt_node->intervals->data[i];
        for(int_t j = interval->lo; j<= interval->hi; j++) {
            fmd_imt_interval_t* single;
            fmd_imt_interval_init(&single,j,j);
            fmd_vector_append(tmp_lists[vertex_last_char_enc[j]], single);
        }
    }
    for(int_t i = 0; i < alphabet_size; i++) {
        fmd_imt_init_compact_intervals(tmp_lists[i], &o->interval_buckets[i]);
        fmd_vector_free(tmp_lists[i]);
    }
    free(tmp_lists);
    o->left = fmd_oimt_init_helper((fmd_imt_node_t*)imt_node->left, vertex_last_char_enc, alphabet, alphabet_size);
    o->right = fmd_oimt_init_helper((fmd_imt_node_t*)imt_node->right, vertex_last_char_enc, alphabet, alphabet_size);
    return o;
}

void fmd_oimt_query(fmd_oimt_t *oimt, int_t start, int_t end, char_t c, int_t no_max_intervals, fmd_vector_t **intervals) {
    fmd_vector_t *merge_list;
    fmd_vector_init(&merge_list, FMD_VECTOR_INIT_SIZE, &prm_fstruct);

    int_t no_cur_intervals = 0;
    fmd_oimt_query_helper(oimt->root, start, end, oimt->c2e[c], no_max_intervals, &no_cur_intervals, merge_list);

    fmd_vector_t *merged = fmd_imt_multiway_merge_intervals(merge_list, no_max_intervals);
    fmd_vector_free(merge_list);
    *intervals = merged;
}

void fmd_oimt_query_helper(fmd_oimt_node_t *node, int_t lo, int_t hi, int_t enc, int_t no_max_intervals, int_t *no_cur_intervals, fmd_vector_t *merge_list) {
    if(no_max_intervals != -1 && *no_cur_intervals >= no_max_intervals) {
        return;
    }
    if(lo == node->lo && hi == node->hi) {
        fmd_vector_append(merge_list, node->interval_buckets[enc]);
        *no_cur_intervals += node->interval_buckets[enc]->size;
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            fmd_oimt_query_helper((fmd_oimt_node_t*)node->left, lo, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        } else if (lo > split) {
            fmd_oimt_query_helper((fmd_oimt_node_t*)node->right, lo, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        } else {
            fmd_oimt_query_helper((fmd_oimt_node_t*)node->left, lo, split, enc, no_max_intervals, no_cur_intervals, merge_list);
            fmd_oimt_query_helper((fmd_oimt_node_t*)node->right, split + 1, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        }
    }
}

void fmd_oimt_free(fmd_oimt_t *oimt) {
    if(!oimt) return;
    fmd_oimt_free_helper(oimt->root, oimt->alphabet_size);
    free(oimt->alphabet);
    free(oimt);
}

void fmd_oimt_free_helper(fmd_oimt_node_t *oimt_node, int_t alphabet_size) {
    if(!oimt_node) return;
    for(int_t i = 0; i < alphabet_size; i++) {
        fmd_vector_free(oimt_node->interval_buckets[i]);
    }
    free(oimt_node->interval_buckets);
    fmd_oimt_node_t *left = oimt_node->left;
    fmd_oimt_node_t *right = oimt_node->right;
    free(oimt_node);
    fmd_oimt_free_helper(left,  alphabet_size);
    fmd_oimt_free_helper(right, alphabet_size);
}

fmd_oimt_t* fmd_oimt_copy(fmd_oimt_t *o);
uint_t fmd_oimt_hash(fmd_oimt_t *o);
void fmd_oimt_comp(fmd_oimt_t *o1, fmd_oimt_t *o2);