#include "fmd_match_chain.h"

void fmd_fmd_match_chain_init(fmd_match_chain_t **chain) {
    *chain = (fmd_match_chain_t*)calloc(1, sizeof(fmd_match_chain_t));
    (*chain)->dummy = (fmd_fmd_match_node_t*)calloc(1, sizeof(fmd_fmd_match_node_t));
    (*chain)->dummy->next = (*chain)->dummy;
    (*chain)->dummy->prev = (*chain)->dummy;
    (*chain)->head = NULL;
    (*chain)->tail = NULL;
    (*chain)->size = 0;
}

void fmd_fmd_match_chain_free(fmd_match_chain_t *chain) {
    fmd_fmd_match_node_t *current = chain->dummy->next;
    while(current != chain->dummy) {
        fmd_fmd_match_node_t *next = current->next;
        fmd_string_free(current->matching_substring);
        free(current);
        current = next;
    }
    free(chain->dummy);
    free(chain);
}

fmd_match_chain_t *fmd_fmd_match_chain_copy(fmd_match_chain_t *chain) {
    fmd_match_chain_t *copy;
    fmd_fmd_match_chain_init(&copy);

    fmd_fmd_match_node_t *current = chain->dummy->next;
    while(current != chain->dummy) {
        fmd_fmd_match_chain_append(copy, current->matching_substring, current->v_lo, current->v_hi, current->sa_lo, current->sa_hi);
        current = current->next;
    }

    return copy;
}

uint_t fmd_fmd_match_chain_hash(fmd_match_chain_t *chain) {
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;
    fmd_fmd_match_node_t *current = chain->dummy->next;
    while(current != chain->dummy) {
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

int fmd_fmd_match_chain_comp(fmd_match_chain_t *l1, fmd_match_chain_t *l2) {
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

void fmd_fmd_match_chain_append(fmd_match_chain_t *chain, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi) {
    fmd_fmd_match_node_t *new_node = (fmd_fmd_match_node_t*)calloc(1,sizeof(fmd_fmd_match_node_t));
    new_node->matching_substring = match_string;
    new_node->v_lo = v_lo;
    new_node->v_hi = v_hi;
    new_node->sa_lo = sa_lo;
    new_node->sa_hi = sa_hi;
    if (chain->size == 0) {
        chain->head = new_node;
        chain->tail = new_node;
        new_node->prev = chain->dummy;
        new_node->next = chain->dummy;
        chain->dummy->prev = new_node;
        chain->dummy->next = new_node;
    } else {
        new_node->prev = chain->tail;
        new_node->next = chain->dummy;
        chain->tail->next = new_node;
        chain->dummy->prev = new_node;
        chain->tail = new_node;
    }
    chain->size++;
}

void fmd_fmd_match_chain_prepend(fmd_match_chain_t *chain, fmd_string_t *match_string, int_t v_lo, int_t v_hi, int_t sa_lo, int_t sa_hi) {
    fmd_fmd_match_node_t *new_node = (fmd_fmd_match_node_t*)calloc(1,sizeof(fmd_fmd_match_node_t));
    new_node->matching_substring = match_string;
    new_node->v_lo = v_lo;
    new_node->v_hi = v_hi;
    new_node->sa_lo = sa_lo;
    new_node->sa_hi = sa_hi;
    if (chain->size == 0) {
        chain->head = new_node;
        chain->tail = new_node;
        new_node->prev = chain->dummy;
        new_node->next = chain->dummy;
        chain->dummy->prev = new_node;
        chain->dummy->next = new_node;
    } else {
        new_node->next = chain->head;
        new_node->prev = chain->dummy;
        chain->head->prev = new_node;
        chain->dummy->next = new_node;
        chain->head = new_node;
    }
    chain->size++;
}
