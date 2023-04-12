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
    if(i){
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
    return (int)(i1->lo - i2->lo);
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

void fmd_imt_query(fmd_imt_t *i, int_t start, int_t end, fmd_vector_t **intervals) {
    fmd_vector_t *ns = fmd_imt_query_helper(i->root, start, end);
    *intervals = ns;
}

fmd_vector_t *fmd_imt_query_helper(fmd_imt_node_t *node, int_t lo, int_t hi) {
    if(node->lo == node->hi) {
        return node->intervals;
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            return fmd_imt_query_helper((fmd_imt_node_t*)node->left, lo, hi);
        } else if (lo > split) {
            return fmd_imt_query_helper((fmd_imt_node_t*)node->right, lo, hi);
        } else {
            return fmd_imt_merge_intervals(fmd_imt_query_helper((fmd_imt_node_t*)node->left, lo, split),
                                           fmd_imt_query_helper((fmd_imt_node_t*)node->right, split+1, hi));
        }
    }
}

fmd_imt_node_t *fmd_imt_init_helper(int_t lo, int_t hi, fmd_vector_t *kv_interval_pairs) {
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
    return result;
}