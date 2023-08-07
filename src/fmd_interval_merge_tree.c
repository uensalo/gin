#include "fmd_interval_merge_tree.h"

void fmd_imt_interval_init(fmd_imt_interval_t **i, int_t lo, int_t hi) {
    fmd_imt_interval_t *it = calloc(1, sizeof(fmd_imt_interval_t));
    if(!it) {
        *i = NULL;
        return;
    }
    it->lo = lo;
    it->hi = hi;
    *i = it;
}

void fmd_imt_interval_free(fmd_imt_interval_t *i) {
    if(i) {
        free(i);
    }
}

fmd_imt_interval_t *fmd_imt_interval_copy(fmd_imt_interval_t *i) {
    if(!i) return NULL;
    fmd_imt_interval_t *c = calloc(1, sizeof(fmd_imt_interval_t));
    if(!c) return NULL;
    c->lo = i->lo;
    c->hi = i->hi;
    return c;
}

uint_t fmd_imt_interval_hash(fmd_imt_interval_t *i) {
    // cantor pairing
    uint_t sum = i->lo + i->hi;
    return (sum * (sum + 1)) / 2 + i->hi;
}

int fmd_imt_interval_comp(fmd_imt_interval_t *i1, fmd_imt_interval_t *i2) {
    return (int)(i1->lo - i2->lo);
}


void fmd_imt_node_init(fmd_imt_node_t **i, int_t lo, int_t hi, fmd_vector_t *intervals, fmd_imt_node_t *left, fmd_imt_node_t *right) {
    fmd_imt_node_t *n = calloc(1, sizeof(fmd_imt_node_t));
    if(!n) {
        *i = NULL;
        return;
    }
    n->lo = lo;
    n->hi = hi;
    n->intervals = intervals;
    n->left = (struct fmd_imt_node_t*)left;
    n->right = (struct fmd_imt_node_t*)right;
    *i = n;
}

void fmd_imt_node_free(fmd_imt_node_t *i) {
    if(i) {
        fmd_vector_free(i->intervals);
        free(i);
    }
}

fmd_imt_node_t *fmd_imt_node_copy(fmd_imt_node_t *i) {
    fmd_imt_node_t *n = calloc(1, sizeof(fmd_imt_node_t));
    if(!n) {
        return NULL;
    }
    n->lo = i->lo;
    n->hi = i->hi;
    n->intervals = fmd_vector_copy(i->intervals);
    n->left = (struct fmd_imt_node_t*)i->left;
    n->right = (struct fmd_imt_node_t*)i->right;
    return n;
}

uint_t fmd_imt_node_hash(fmd_imt_node_t *i) {
    // cantor pairing
    uint_t sum = i->lo + i->hi;
    return (sum * (sum + 1)) / 2 + i->hi;
}

int fmd_imt_node_comp(fmd_imt_node_t *i1, fmd_imt_node_t *i2) {
    return (int)(i1->lo - i2->lo) == 0 ? fmd_vector_comp(i1->intervals, i2->intervals) : -1;
}

void fmd_imt_init(fmd_imt_t **i, int_t no_keys, fmd_vector_t *kv_interval_pairs) {
    fmd_imt_t *t = calloc(1, sizeof(fmd_imt_t));
    if(!t) {
        *i = NULL;
        return;
    }
    t->root = fmd_imt_init_helper(0, no_keys-1, kv_interval_pairs);
    if(!t->root) {
        free(t);
        *i = NULL;
        return;
    }
    t->no_keys = no_keys;
    *i = t;
}

void fmd_imt_free(fmd_imt_t *i) {
    if(i) {
        fmd_imt_free_helper(i->root);
        free(i);
    }
}

fmd_imt_node_t *fmd_imt_copy_helper(fmd_imt_node_t *i) {
    if (!i) return NULL;
    fmd_imt_node_t *n = calloc(1, sizeof(fmd_imt_node_t));
    if (!n) return NULL;
    n->lo = i->lo;
    n->hi = i->hi;
    n->intervals = fmd_vector_copy(i->intervals);
    n->left = (struct fmd_imt_node_t*)fmd_imt_copy_helper((fmd_imt_node_t*)i->left);
    n->right = (struct fmd_imt_node_t*)fmd_imt_copy_helper((fmd_imt_node_t*)i->right);
    return n;
}

uint_t fmd_imt_hash_helper(fmd_imt_node_t *i) {
    if (!i) return 0;
    uint_t hash_val = fmd_imt_node_hash(i);
    hash_val ^= fmd_imt_hash_helper((fmd_imt_node_t*)i->left) * 31;
    hash_val ^= fmd_imt_hash_helper((fmd_imt_node_t*)i->right) * 37;
    return hash_val;
}

int fmd_imt_comp_helper(fmd_imt_node_t *i1, fmd_imt_node_t *i2) {
    if (!i1 && !i2) return 0;
    if (!i1) return -1;
    if (!i2) return 1;
    int cmp_val = fmd_imt_node_comp(i1, i2);
    if (cmp_val != 0) return cmp_val;
    cmp_val = fmd_imt_comp_helper((fmd_imt_node_t*)i1->left, (fmd_imt_node_t*)i2->left);
    if (cmp_val != 0) return cmp_val;
    return fmd_imt_comp_helper((fmd_imt_node_t*)i1->right, (fmd_imt_node_t*)i2->right);
}

fmd_imt_t *fmd_imt_copy(fmd_imt_t *i) {
    if (!i) return NULL;
    fmd_imt_t *copy = calloc(1, sizeof(fmd_imt_t));
    if (!copy) return NULL;
    copy->no_keys = i->no_keys;
    copy->root = fmd_imt_copy_helper(i->root);
    return copy;
}

uint_t fmd_imt_hash(fmd_imt_t *i) {
    if (!i) return 0;
    return fmd_imt_hash_helper(i->root);
}

int fmd_imt_comp(fmd_imt_t *i1, fmd_imt_t *i2) {
    if (!i1 && !i2) return 0;
    if (!i1) return -1;
    if (!i2) return 1;
    return fmd_imt_comp_helper(i1->root, i2->root);
}

void fmd_imt_query(fmd_imt_t *i, int_t start, int_t end, int_t no_max_intervals, fmd_vector_t **intervals) {
    fmd_vector_t *merge_list;
    fmd_vector_init(&merge_list, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    // gather phase
    int_t no_cur_intervals = 0;
    fmd_imt_query_helper(i->root, start, end, no_max_intervals, &no_cur_intervals, merge_list);

    // k-way merge phase
    fmd_vector_t *merged = fmd_imt_multiway_merge_intervals(merge_list, no_max_intervals);
    fmd_vector_free(merge_list);
    *intervals = merged;
}

void fmd_imt_query_helper(fmd_imt_node_t *node, int_t lo, int_t hi, int_t no_max_intervals, int_t *no_cur_intervals, fmd_vector_t *merge_list) {
    if(no_max_intervals != -1 && *no_cur_intervals >= no_max_intervals) {
        return;
    }
    if(node->lo == node->hi) { // lo == node->lo && hi == node->hi
        fmd_vector_append(merge_list, node->intervals); // important to return a copy for memory management purposes
        *no_cur_intervals += node->intervals->size;
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            fmd_imt_query_helper((fmd_imt_node_t*)node->left, lo, hi, no_max_intervals, no_cur_intervals, merge_list);
        } else if (lo > split) {
            fmd_imt_query_helper((fmd_imt_node_t*)node->right, lo, hi, no_max_intervals, no_cur_intervals, merge_list);
        } else {
            fmd_imt_query_helper((fmd_imt_node_t *) node->left, lo, split, no_max_intervals, no_cur_intervals, merge_list);
            fmd_imt_query_helper((fmd_imt_node_t *) node->right, split + 1, hi, no_max_intervals, no_cur_intervals, merge_list);
        }
    }
}

void fmd_imt_query_legacy(fmd_imt_t *i, int_t start, int_t end, fmd_vector_t **intervals) {
    fmd_vector_t *ns = fmd_imt_query_helper_legacy(i->root, start, end);
    *intervals = ns;
}

fmd_vector_t *fmd_imt_query_helper_legacy(fmd_imt_node_t *node, int_t lo, int_t hi) {
    if(node->lo == node->hi) {
        return fmd_vector_copy(node->intervals); // important to return a copy for memory management purposes
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            return fmd_imt_query_helper_legacy((fmd_imt_node_t*)node->left, lo, hi);
        } else if (lo > split) {
            return fmd_imt_query_helper_legacy((fmd_imt_node_t*)node->right, lo, hi);
        } else {
            fmd_vector_t *i1, *i2, *merged;
            #pragma omp task default(none) shared(node, lo, split, hi, i1)
            {
                i1 = fmd_imt_query_helper_legacy((fmd_imt_node_t *) node->left, lo, split);
            }
            #pragma omp task default(none) shared(node, lo, split, hi, i2)
            {
                i2 = fmd_imt_query_helper_legacy((fmd_imt_node_t *) node->right, split + 1, hi);
            }
            #pragma omp taskwait
            merged = fmd_imt_merge_intervals(i1, i2);
            return merged;
        }
    }
}

void fmd_imt_init_compact_intervals(fmd_vector_t *intervals, fmd_vector_t **compacted) {
    fmd_vector_sort(intervals);
    fmd_vector_init(compacted, intervals->size, &fmd_fstruct_imt_interval);
    if(intervals->size) {
        fmd_imt_interval_t *cur_interval = fmd_imt_interval_copy(intervals->data[0]);
        fmd_vector_append(*compacted, cur_interval);
        for (int_t i = 1; i < intervals->size; i++) {
            fmd_imt_interval_t *next_interval = intervals->data[i];
            if(cur_interval->hi + 1 >= next_interval->lo) {
                cur_interval->hi = next_interval->hi > cur_interval->hi ? next_interval->hi : cur_interval->hi;
            } else {
                cur_interval = fmd_imt_interval_copy(next_interval);
                fmd_vector_append(*compacted, cur_interval);
            }
        }
    }
    fmd_vector_fit(*compacted);
}

fmd_imt_node_t *fmd_imt_init_helper(int_t lo, int_t hi, fmd_vector_t *kv_interval_pairs) {
    fmd_imt_node_t *n = calloc(1, sizeof(fmd_imt_node_t));
    n->lo = lo;
    n->hi = hi;

    // compute the pre merged list
    fmd_vector_init(&n->intervals, hi - lo + 1, &fmd_fstruct_imt_interval);
    for(int_t i = lo; i <= hi; i++) {
        fmd_vector_t *value_intervals = (fmd_vector_t*)kv_interval_pairs->data[i];
        for(int_t j = 0; j < value_intervals->size; j++) {
            fmd_imt_interval_t *interval = fmd_imt_interval_copy(value_intervals->data[j]);
            fmd_vector_append(n->intervals, interval);
        }
    }
    //fmd_vector_sort(n->intervals);
    fmd_vector_t *compacted;
    fmd_imt_init_compact_intervals(n->intervals, &compacted);
    fmd_vector_free(n->intervals);
    n->intervals = compacted;

    if(lo != hi) {
        int_t split = (lo + hi) / 2;
        n->left = (struct fmd_imt_node_t*)fmd_imt_init_helper(lo, split, kv_interval_pairs);
        n->right = (struct fmd_imt_node_t*)fmd_imt_init_helper(split+1, hi, kv_interval_pairs);
    }
    return n;


    /*
    fmd_imt_node_t *n = calloc(1, sizeof(fmd_imt_node_t));
    n->lo = lo;
    n->hi = hi;
    if(lo == hi) {
        n->intervals = kv_interval_pairs->data[lo];
    } else {
        int_t split = (lo + hi) / 2;
        n->left = (struct fmd_imt_node_t*)fmd_imt_init_helper(lo, split, kv_interval_pairs);
        n->right = (struct fmd_imt_node_t*)fmd_imt_init_helper(split+1, hi, kv_interval_pairs);
    }
    return n;
    */
}

void fmd_imt_free_helper(fmd_imt_node_t *node) {
    if(node) {
        fmd_imt_free_helper((fmd_imt_node_t*)node->left);
        fmd_imt_free_helper((fmd_imt_node_t*)node->right);
        fmd_imt_node_free(node);
    }
}

fmd_vector_t *fmd_imt_merge_intervals(fmd_vector_t *i1, fmd_vector_t *i2) {
    fmd_vector_t *result;
    fmd_vector_init(&result, i1->size + i2->size, &fmd_fstruct_imt_interval);

    int_t index1 = 0, index2 = 0;
    fmd_imt_interval_t *current = NULL, *next = NULL;

    while (index1 < i1->size || index2 < i2->size) {
        if (index1 < i1->size && (index2 == i2->size || fmd_imt_interval_comp(i1->data[index1], i2->data[index2]) < 0)) {
            next = i1->data[index1++];
        } else {
            next = i2->data[index2++];
        }
        if (!current) {
            current = fmd_imt_interval_copy(next);
        } else if (current->hi + 1 >= next->lo) {
            if (next->hi > current->hi) {
                current->hi = next->hi;
            }
        } else {
            fmd_vector_append(result, current);
            current = fmd_imt_interval_copy(next);
        }
    }
    if (current) {
        fmd_vector_append(result, current);
    }
    fmd_vector_free(i1);
    fmd_vector_free(i2);
    return result;
}

fmd_vector_t *fmd_imt_multiway_merge_intervals(fmd_vector_t *list_of_intervals, int_t no_max_intervals) {
    fmd_vector_t *merged;
    fmd_vector_init(&merged, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_imt_interval);
    if(list_of_intervals->size == 1) {
        fmd_vector_t *intervals = list_of_intervals->data[0];
        int_t bound = no_max_intervals == -1 ? intervals->size : MIN2(no_max_intervals, intervals->size);
        for(int_t i = 0; i < bound; i++) {
            fmd_imt_interval_t *interval = fmd_imt_interval_copy(intervals->data[i]);
            fmd_vector_append(merged,interval);
        }
        return merged;
    } else {
        fmd_min_heap_t *pq; // key: intervals, value: index of interval list
        fmd_min_heap_init(&pq, list_of_intervals->size, &fmd_fstruct_imt_interval, &prm_fstruct);

        int_t *ptrs = calloc(list_of_intervals->size, sizeof(int_t));
        for (int_t i = 0; i < list_of_intervals->size; i++) {
            fmd_vector_t *interval_list = list_of_intervals->data[i];
            if (interval_list->size) {
                fmd_min_heap_push(pq, interval_list->data[ptrs[i]++], (void *) i);
            }
        }

        fmd_imt_interval_t *cur = NULL, *next = NULL;
        while (pq->size && !(no_max_intervals != -1 && merged->size >= no_max_intervals)) {
            int_t idx;
            fmd_min_heap_pop(pq, &next, &idx);
            if (!cur) {
                cur = fmd_imt_interval_copy(next);
            } else if (cur->hi + 1 >= next->lo) {
                if (next->hi > cur->hi) {
                    cur->hi = next->hi;
                }
            } else {
                fmd_vector_append(merged, cur);
                cur = fmd_imt_interval_copy(next);
            }
            fmd_vector_t *interval_list = list_of_intervals->data[idx];
            if (ptrs[idx] < interval_list->size) {
                fmd_min_heap_push(pq, interval_list->data[ptrs[idx]++], (void *) idx);
            }
        }
        if (cur) {
            fmd_vector_append(merged, cur);
        }
        free(ptrs);
        fmd_min_heap_free(pq);
        return merged;
    }
}