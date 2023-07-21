#include "fmd_match_list.h"

void fmd_fmd_match_list_init(fmd_fmd_match_list_t **list) {
    *list = (fmd_fmd_match_list_t*)calloc(1, sizeof(fmd_fmd_match_list_t));
    (*list)->dummy = (fmd_fmd_match_node_t*)calloc(1, sizeof(fmd_fmd_match_node_t));
    (*list)->dummy->next = (*list)->dummy;
    (*list)->dummy->prev = (*list)->dummy;
    (*list)->head = NULL;
    (*list)->tail = NULL;
    (*list)->size = 0;
}

void fmd_fmd_match_list_free(fmd_fmd_match_list_t *list) {
    fmd_fmd_match_node_t *current = list->dummy->next;
    while(current != list->dummy) {
        fmd_fmd_match_node_t *next = current->next;
        fmd_string_free(current->matching_substring);
        free(current);
        current = next;
    }
    free(list->dummy);
    free(list);
}

fmd_fmd_match_list_t *fmd_fmd_match_list_copy(fmd_fmd_match_list_t *list) {
    fmd_fmd_match_list_t *copy;
    fmd_fmd_match_list_init(&copy);

    fmd_fmd_match_node_t *current = list->dummy->next;
    while(current != list->dummy) {
        fmd_fmd_match_list_append(copy, current->matching_substring, current->v_lo, current->v_hi, current->sa_lo, current->sa_hi);
        current = current->next;
    }

    return copy;
}

uint_t fmd_fmd_match_list_hash(fmd_fmd_match_list_t *list) {
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;
    fmd_fmd_match_node_t *current = list->dummy->next;
    while(current != list->dummy) {
        uint_t hash_item = fmd_string_hash(current->matching_substring)
                ^ prm_hash_f((void*)current->v_lo)
                ^ prm_hash_f((void*)current->v_hi)
                ^ prm_hash_f((void*)current->sa_lo)
                ^ prm_hash_f((void*)current->sa_hi);
        hash ^= hash_item;
        hash *= prime;
        current = current->next;
    }
    return hash;
}

int fmd_fmd_match_list_comp(fmd_fmd_match_list_t *l1, fmd_fmd_match_list_t *l2) {
    fmd_fmd_match_node_t *n1 = l1->dummy->next;
    fmd_fmd_match_node_t *n2 = l2->dummy->next;
    while(n1 != l1->dummy && n2 != l2->dummy) {
        int cmp = fmd_string_comp(n1->matching_substring, n2->matching_substring);
        if(cmp != 0) {
            return cmp;
        }
        n1 = n1->next;
        n2 = n2->next;
    }
    return (n1 != l1->dummy) - (n2 != l2->dummy);
}

void fmd_fmd_match_list_append(fmd_fmd_match_list_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi) {
    fmd_fmd_match_node_t *new_node = (fmd_fmd_match_node_t*)calloc(1,sizeof(fmd_fmd_match_node_t));
    new_node->matching_substring = match_string;
    new_node->v_lo = v_lo;
    new_node->v_hi = v_hi;
    new_node->sa_lo = sa_lo;
    new_node->sa_hi = sa_hi;

    if (list->size == 0) {
        list->head = new_node;
        list->tail = new_node;
        new_node->prev = list->dummy;
        new_node->next = list->dummy;
        list->dummy->prev = new_node;
        list->dummy->next = new_node;
    } else {
        new_node->prev = list->tail;
        new_node->next = list->dummy;
        list->tail->next = new_node;
        list->dummy->prev = new_node;
        list->tail = new_node;
    }
    list->size++;
}

void fmd_fmd_match_list_prepend(fmd_fmd_match_list_t *list, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi) {
    fmd_fmd_match_node_t *new_node = (fmd_fmd_match_node_t*)calloc(1,sizeof(fmd_fmd_match_node_t));
    new_node->matching_substring = match_string;
    new_node->v_lo = v_lo;
    new_node->v_hi = v_hi;
    new_node->sa_lo = sa_lo;
    new_node->sa_hi = sa_hi;

    if (list->size == 0) {
        list->head = new_node;
        list->tail = new_node;
        new_node->prev = list->dummy;
        new_node->next = list->dummy;
        list->dummy->prev = new_node;
        list->dummy->next = new_node;
    } else {
        new_node->next = list->head;
        new_node->prev = list->dummy;
        list->head->prev = new_node;
        list->dummy->next = new_node;
        list->head = new_node;
    }
    list->size++;
}
