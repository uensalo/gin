#include "gin_oimt.h"

void gin_oimt_init(gin_imt_t *imt, int_t *vertex_last_char_enc, int_t *alphabet, int_t alphabet_size, gin_oimt_t **oimt) {
    gin_oimt_t *o = calloc(1, sizeof(gin_oimt_t));
    o->root = gin_oimt_init_helper(imt->root, vertex_last_char_enc, alphabet, alphabet_size);
    o->no_keys = imt->no_keys;
    o->alphabet = calloc(alphabet_size, sizeof(int_t));
    memcpy(o->alphabet, alphabet, alphabet_size * sizeof(int_t));
    o->alphabet_size = alphabet_size;

    int_t t = 0;
    o->c2e = calloc(GIN_MAX_ALPHABET_SIZE, sizeof(int_t));
    for(int_t i = 0; i < alphabet_size; i++) {
        o->c2e[alphabet[i]] = t;
        t++;
    }
    *oimt = o;
}

gin_oimt_node_t *gin_oimt_init_helper(gin_imt_node_t *imt_node,
                                      int_t *vertex_last_char_enc,
                                      int_t *alphabet,
                                      int_t alphabet_size) {
    if(!imt_node) return NULL;
    gin_oimt_node_t *o = calloc(1, sizeof(gin_oimt_node_t));
    o->lo = imt_node->lo;
    o->hi = imt_node->hi;
    // partition the alphabet
    o->interval_buckets = calloc(alphabet_size, sizeof(gin_vector_t*));
    gin_vector_t **tmp_lists = calloc(alphabet_size, sizeof(gin_vector_t*));
    for(int_t i = 0; i < alphabet_size; i++) {
        gin_vector_init(&tmp_lists[i], GIN_VECTOR_INIT_SIZE, &gin_fstruct_imt_interval);
    }

    for(int_t i = 0; i < imt_node->intervals->size; i++) {
        gin_imt_interval_t *interval = imt_node->intervals->data[i];
        for(int_t j = interval->lo; j<= interval->hi; j++) {
            gin_imt_interval_t* single;
            gin_imt_interval_init(&single,j,j);
            gin_vector_append(tmp_lists[vertex_last_char_enc[j]], single);
        }
    }
    for(int_t i = 0; i < alphabet_size; i++) {
        gin_imt_init_compact_intervals(tmp_lists[i], &o->interval_buckets[i]);
        gin_vector_free(tmp_lists[i]);
    }
    free(tmp_lists);
    o->left = gin_oimt_init_helper((gin_imt_node_t*)imt_node->left, vertex_last_char_enc, alphabet, alphabet_size);
    o->right = gin_oimt_init_helper((gin_imt_node_t*)imt_node->right, vertex_last_char_enc, alphabet, alphabet_size);
    return o;
}

void gin_oimt_query(gin_oimt_t *oimt, int_t start, int_t end, char_t c, int_t no_max_intervals, gin_vector_t **intervals) {
    gin_vector_t *merge_list;
    gin_vector_init(&merge_list, GIN_VECTOR_INIT_SIZE, &prm_fstruct);

    int_t no_cur_intervals = 0;
    gin_oimt_query_helper(oimt->root, start, end, oimt->c2e[c], no_max_intervals, &no_cur_intervals, merge_list);

    gin_vector_t *merged = gin_imt_multiway_merge_intervals(merge_list, no_max_intervals);
    gin_vector_free(merge_list);
    *intervals = merged;
}

void gin_oimt_query_helper(gin_oimt_node_t *node, int_t lo, int_t hi, int_t enc, int_t no_max_intervals, int_t *no_cur_intervals, gin_vector_t *merge_list) {
    if(no_max_intervals != -1 && *no_cur_intervals >= no_max_intervals) {
        return;
    }
    if(lo == node->lo && hi == node->hi) {
        gin_vector_append(merge_list, node->interval_buckets[enc]);
        *no_cur_intervals += node->interval_buckets[enc]->size;
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            gin_oimt_query_helper((gin_oimt_node_t*)node->left, lo, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        } else if (lo > split) {
            gin_oimt_query_helper((gin_oimt_node_t*)node->right, lo, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        } else {
            gin_oimt_query_helper((gin_oimt_node_t*)node->left, lo, split, enc, no_max_intervals, no_cur_intervals, merge_list);
            gin_oimt_query_helper((gin_oimt_node_t*)node->right, split + 1, hi, enc, no_max_intervals, no_cur_intervals, merge_list);
        }
    }
}

void gin_oimt_free(gin_oimt_t *oimt) {
    if(!oimt) return;
    gin_oimt_free_helper(oimt->root, oimt->alphabet_size);
    free(oimt->alphabet);
    free(oimt);
}

void gin_oimt_free_helper(gin_oimt_node_t *oimt_node, int_t alphabet_size) {
    if(!oimt_node) return;
    for(int_t i = 0; i < alphabet_size; i++) {
        gin_vector_free(oimt_node->interval_buckets[i]);
    }
    free(oimt_node->interval_buckets);
    gin_oimt_node_t *left = oimt_node->left;
    gin_oimt_node_t *right = oimt_node->right;
    free(oimt_node);
    gin_oimt_free_helper(left,  alphabet_size);
    gin_oimt_free_helper(right, alphabet_size);
}

gin_oimt_t* gin_oimt_copy(gin_oimt_t *o);
uint_t gin_oimt_hash(gin_oimt_t *o);
void gin_oimt_comp(gin_oimt_t *o1, gin_oimt_t *o2);