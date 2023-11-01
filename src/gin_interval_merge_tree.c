/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_interval_merge_tree.c is part of gin
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
#include "gin_interval_merge_tree.h"

void gin_imt_interval_init(gin_imt_interval_t **i, int_t lo, int_t hi) {
    gin_imt_interval_t *it = calloc(1, sizeof(gin_imt_interval_t));
    if(!it) {
        *i = NULL;
        return;
    }
    it->lo = lo;
    it->hi = hi;
    *i = it;
}

void gin_imt_interval_free(gin_imt_interval_t *i) {
    if(i) {
        free(i);
    }
}

gin_imt_interval_t *gin_imt_interval_copy(gin_imt_interval_t *i) {
    if(!i) return NULL;
    gin_imt_interval_t *c = calloc(1, sizeof(gin_imt_interval_t));
    if(!c) return NULL;
    c->lo = i->lo;
    c->hi = i->hi;
    return c;
}

uint_t gin_imt_interval_hash(gin_imt_interval_t *i) {
    // cantor pairing
    uint_t sum = i->lo + i->hi;
    return (sum * (sum + 1)) / 2 + i->hi;
}

int gin_imt_interval_comp(gin_imt_interval_t *i1, gin_imt_interval_t *i2) {
    return (int)(i1->lo - i2->lo);
}


void gin_imt_node_init(gin_imt_node_t **i, int_t lo, int_t hi, gin_vector_t *intervals, gin_imt_node_t *left, gin_imt_node_t *right) {
    gin_imt_node_t *n = calloc(1, sizeof(gin_imt_node_t));
    if(!n) {
        *i = NULL;
        return;
    }
    n->lo = lo;
    n->hi = hi;
    n->intervals = intervals;
    n->left = (struct gin_imt_node_t*)left;
    n->right = (struct gin_imt_node_t*)right;
    *i = n;
}

void gin_imt_node_free(gin_imt_node_t *i) {
    if(i) {
        gin_vector_free(i->intervals);
        free(i);
    }
}

gin_imt_node_t *gin_imt_node_copy(gin_imt_node_t *i) {
    gin_imt_node_t *n = calloc(1, sizeof(gin_imt_node_t));
    if(!n) {
        return NULL;
    }
    n->lo = i->lo;
    n->hi = i->hi;
    n->intervals = gin_vector_copy(i->intervals);
    n->left = (struct gin_imt_node_t*)i->left;
    n->right = (struct gin_imt_node_t*)i->right;
    return n;
}

uint_t gin_imt_node_hash(gin_imt_node_t *i) {
    // cantor pairing
    uint_t sum = i->lo + i->hi;
    return (sum * (sum + 1)) / 2 + i->hi;
}

int gin_imt_node_comp(gin_imt_node_t *i1, gin_imt_node_t *i2) {
    return (int)(i1->lo - i2->lo) == 0 ? gin_vector_comp(i1->intervals, i2->intervals) : -1;
}

void gin_imt_init(gin_imt_t **i, int_t no_keys, gin_vector_t *kv_interval_pairs) {
    gin_imt_t *t = calloc(1, sizeof(gin_imt_t));
    if(!t) {
        *i = NULL;
        return;
    }
    t->root = gin_imt_init_helper(0, no_keys-1, kv_interval_pairs);
    if(!t->root) {
        free(t);
        *i = NULL;
        return;
    }
    t->no_keys = no_keys;
    *i = t;
}

void gin_imt_free(gin_imt_t *i) {
    if(i) {
        gin_imt_free_helper(i->root);
        free(i);
    }
}

gin_imt_node_t *gin_imt_copy_helper(gin_imt_node_t *i) {
    if (!i) return NULL;
    gin_imt_node_t *n = calloc(1, sizeof(gin_imt_node_t));
    if (!n) return NULL;
    n->lo = i->lo;
    n->hi = i->hi;
    n->intervals = gin_vector_copy(i->intervals);
    n->left = (struct gin_imt_node_t*)gin_imt_copy_helper((gin_imt_node_t*)i->left);
    n->right = (struct gin_imt_node_t*)gin_imt_copy_helper((gin_imt_node_t*)i->right);
    return n;
}

uint_t gin_imt_hash_helper(gin_imt_node_t *i) {
    if (!i) return 0;
    uint_t hash_val = gin_imt_node_hash(i);
    hash_val ^= gin_imt_hash_helper((gin_imt_node_t*)i->left) * 31;
    hash_val ^= gin_imt_hash_helper((gin_imt_node_t*)i->right) * 37;
    return hash_val;
}

int gin_imt_comp_helper(gin_imt_node_t *i1, gin_imt_node_t *i2) {
    if (!i1 && !i2) return 0;
    if (!i1) return -1;
    if (!i2) return 1;
    int cmp_val = gin_imt_node_comp(i1, i2);
    if (cmp_val != 0) return cmp_val;
    cmp_val = gin_imt_comp_helper((gin_imt_node_t*)i1->left, (gin_imt_node_t*)i2->left);
    if (cmp_val != 0) return cmp_val;
    return gin_imt_comp_helper((gin_imt_node_t*)i1->right, (gin_imt_node_t*)i2->right);
}

gin_imt_t *gin_imt_copy(gin_imt_t *i) {
    if (!i) return NULL;
    gin_imt_t *copy = calloc(1, sizeof(gin_imt_t));
    if (!copy) return NULL;
    copy->no_keys = i->no_keys;
    copy->root = gin_imt_copy_helper(i->root);
    return copy;
}

uint_t gin_imt_hash(gin_imt_t *i) {
    if (!i) return 0;
    return gin_imt_hash_helper(i->root);
}

int gin_imt_comp(gin_imt_t *i1, gin_imt_t *i2) {
    if (!i1 && !i2) return 0;
    if (!i1) return -1;
    if (!i2) return 1;
    return gin_imt_comp_helper(i1->root, i2->root);
}

void gin_imt_query(gin_imt_t *i, int_t start, int_t end, int_t no_max_intervals, gin_vector_t **intervals) {
    gin_vector_t *merge_list;
    gin_vector_init(&merge_list, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    // gather phase
    int_t no_cur_intervals = 0;
    gin_imt_query_helper(i->root, start, end, no_max_intervals, &no_cur_intervals, merge_list);

    // k-way merge phase
    gin_vector_t *merged = gin_imt_multiway_merge_intervals(merge_list, no_max_intervals);
    gin_vector_free(merge_list);
    *intervals = merged;
}

void gin_imt_query_helper(gin_imt_node_t *node, int_t lo, int_t hi, int_t no_max_intervals, int_t *no_cur_intervals, gin_vector_t *merge_list) {
    if(no_max_intervals != -1 && *no_cur_intervals >= no_max_intervals) {
        return;
    }
    if(lo == node->lo && hi == node->hi) {
        gin_vector_append(merge_list, node->intervals); // important to return a copy for memory management purposes
        *no_cur_intervals += node->intervals->size;
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            gin_imt_query_helper((gin_imt_node_t*)node->left, lo, hi, no_max_intervals, no_cur_intervals, merge_list);
        } else if (lo > split) {
            gin_imt_query_helper((gin_imt_node_t*)node->right, lo, hi, no_max_intervals, no_cur_intervals, merge_list);
        } else {
            gin_imt_query_helper((gin_imt_node_t *) node->left, lo, split, no_max_intervals, no_cur_intervals, merge_list);
            gin_imt_query_helper((gin_imt_node_t *) node->right, split + 1, hi, no_max_intervals, no_cur_intervals, merge_list);
        }
    }
}

void gin_imt_query_legacy(gin_imt_t *i, int_t start, int_t end, gin_vector_t **intervals) {
    gin_vector_t *ns = gin_imt_query_helper_legacy(i->root, start, end);
    *intervals = ns;
}

gin_vector_t *gin_imt_query_helper_legacy(gin_imt_node_t *node, int_t lo, int_t hi) {
    if(node->lo == node->hi) {
        return gin_vector_copy(node->intervals); // important to return a copy for memory management purposes
    } else {
        int_t split = (node->lo + node->hi) / 2;
        if(hi <= split) {
            return gin_imt_query_helper_legacy((gin_imt_node_t*)node->left, lo, hi);
        } else if (lo > split) {
            return gin_imt_query_helper_legacy((gin_imt_node_t*)node->right, lo, hi);
        } else {
            gin_vector_t *i1, *i2, *merged;
            #pragma omp task default(none) shared(node, lo, split, hi, i1)
            {
                i1 = gin_imt_query_helper_legacy((gin_imt_node_t *) node->left, lo, split);
            }
            #pragma omp task default(none) shared(node, lo, split, hi, i2)
            {
                i2 = gin_imt_query_helper_legacy((gin_imt_node_t *) node->right, split + 1, hi);
            }
            #pragma omp taskwait
            merged = gin_imt_merge_intervals(i1, i2);
            return merged;
        }
    }
}

void gin_imt_init_compact_intervals(gin_vector_t *intervals, gin_vector_t **compacted) {
    gin_vector_sort(intervals);
    gin_vector_init(compacted, intervals->size, &gin_fstruct_imt_interval);
    if(intervals->size) {
        gin_imt_interval_t *cur_interval = gin_imt_interval_copy(intervals->data[0]);
        gin_vector_append(*compacted, cur_interval);
        for (int_t i = 1; i < intervals->size; i++) {
            gin_imt_interval_t *next_interval = intervals->data[i];
            if(cur_interval->hi + 1 >= next_interval->lo) {
                cur_interval->hi = next_interval->hi > cur_interval->hi ? next_interval->hi : cur_interval->hi;
            } else {
                cur_interval = gin_imt_interval_copy(next_interval);
                gin_vector_append(*compacted, cur_interval);
            }
        }
    }
    gin_vector_fit(*compacted);
}

gin_imt_node_t *gin_imt_init_helper(int_t lo, int_t hi, gin_vector_t *kv_interval_pairs) {
    gin_imt_node_t *n = calloc(1, sizeof(gin_imt_node_t));
    n->lo = lo;
    n->hi = hi;

    // compute the pre merged list
    gin_vector_init(&n->intervals, hi - lo + 1, &gin_fstruct_imt_interval);
    for(int_t i = lo; i <= hi; i++) {
        gin_vector_t *value_intervals = (gin_vector_t*)kv_interval_pairs->data[i];
        for(int_t j = 0; j < value_intervals->size; j++) {
            gin_imt_interval_t *interval = gin_imt_interval_copy(value_intervals->data[j]);
            gin_vector_append(n->intervals, interval);
        }
    }
    //gin_vector_sort(n->intervals);
    gin_vector_t *compacted;
    gin_imt_init_compact_intervals(n->intervals, &compacted);
    gin_vector_free(n->intervals);
    n->intervals = compacted;

    if(lo != hi) {
        int_t split = (lo + hi) / 2;
        n->left = (struct gin_imt_node_t*)gin_imt_init_helper(lo, split, kv_interval_pairs);
        n->right = (struct gin_imt_node_t*)gin_imt_init_helper(split+1, hi, kv_interval_pairs);
    }
    return n;


    /*
    gin_imt_node_t *n = calloc(1, sizeof(gin_imt_node_t));
    n->lo = lo;
    n->hi = hi;
    if(lo == hi) {
        n->intervals = kv_interval_pairs->data[lo];
    } else {
        int_t split = (lo + hi) / 2;
        n->left = (struct gin_imt_node_t*)gin_imt_init_helper(lo, split, kv_interval_pairs);
        n->right = (struct gin_imt_node_t*)gin_imt_init_helper(split+1, hi, kv_interval_pairs);
    }
    return n;
    */
}

void gin_imt_free_helper(gin_imt_node_t *node) {
    if(node) {
        gin_imt_free_helper((gin_imt_node_t*)node->left);
        gin_imt_free_helper((gin_imt_node_t*)node->right);
        gin_imt_node_free(node);
    }
}

gin_vector_t *gin_imt_merge_intervals(gin_vector_t *i1, gin_vector_t *i2) {
    gin_vector_t *result;
    gin_vector_init(&result, i1->size + i2->size, &gin_fstruct_imt_interval);

    int_t index1 = 0, index2 = 0;
    gin_imt_interval_t *current = NULL, *next = NULL;

    while (index1 < i1->size || index2 < i2->size) {
        if (index1 < i1->size && (index2 == i2->size || gin_imt_interval_comp(i1->data[index1], i2->data[index2]) < 0)) {
            next = i1->data[index1++];
        } else {
            next = i2->data[index2++];
        }
        if (!current) {
            current = gin_imt_interval_copy(next);
        } else if (current->hi + 1 >= next->lo) {
            if (next->hi > current->hi) {
                current->hi = next->hi;
            }
        } else {
            gin_vector_append(result, current);
            current = gin_imt_interval_copy(next);
        }
    }
    if (current) {
        gin_vector_append(result, current);
    }
    gin_vector_free(i1);
    gin_vector_free(i2);
    return result;
}

gin_vector_t *gin_imt_multiway_merge_intervals(gin_vector_t *list_of_intervals, int_t no_max_intervals) {
    gin_vector_t *merged;
    gin_vector_init(&merged, GIN_VECTOR_INIT_SIZE, &gin_fstruct_imt_interval);
    if(list_of_intervals->size == 1) {
        gin_vector_t *intervals = list_of_intervals->data[0];
        int_t bound = no_max_intervals == -1 ? intervals->size : MIN2(no_max_intervals, intervals->size);
        for(int_t i = 0; i < bound; i++) {
            gin_imt_interval_t *interval = gin_imt_interval_copy(intervals->data[i]);
            gin_vector_append(merged,interval);
        }
        return merged;
    } else {
        gin_min_heap_t *pq; // key: intervals, value: index of interval list
        gin_min_heap_init(&pq, list_of_intervals->size, &gin_fstruct_imt_interval, &prm_fstruct);

        int_t *ptrs = calloc(list_of_intervals->size, sizeof(int_t));
        for (int_t i = 0; i < list_of_intervals->size; i++) {
            gin_vector_t *interval_list = list_of_intervals->data[i];
            if (interval_list->size) {
                gin_min_heap_push(pq, interval_list->data[ptrs[i]++], (void *) i);
            }
        }

        gin_imt_interval_t *cur = NULL, *next = NULL;
        while (pq->size && !(no_max_intervals != -1 && merged->size >= no_max_intervals)) {
            int_t idx;
            gin_min_heap_pop(pq, &next, &idx);
            if (!cur) {
                cur = gin_imt_interval_copy(next);
            } else if (cur->hi + 1 >= next->lo) {
                if (next->hi > cur->hi) {
                    cur->hi = next->hi;
                }
            } else {
                gin_vector_append(merged, cur);
                cur = gin_imt_interval_copy(next);
            }
            gin_vector_t *interval_list = list_of_intervals->data[idx];
            if (ptrs[idx] < interval_list->size) {
                gin_min_heap_push(pq, interval_list->data[ptrs[idx]++], (void *) idx);
            }
        }
        if (cur) {
            gin_vector_append(merged, cur);
        }
        free(ptrs);
        gin_min_heap_free(pq);
        return merged;
    }
}