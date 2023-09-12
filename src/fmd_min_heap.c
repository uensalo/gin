/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_min_heap.c is part of fmd
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
#include "fmd_min_heap.h"

void fmd_min_heap_init(fmd_min_heap_t **heap, int_t capacity, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    *heap = (fmd_min_heap_t *)malloc(sizeof(fmd_min_heap_t));
    (*heap)->data = (fmd_min_heap_kv_t *)calloc(capacity, sizeof(fmd_min_heap_kv_t));
    (*heap)->capacity = capacity;
    (*heap)->size = 0;
    (*heap)->key_f = key_f;
    (*heap)->val_f = val_f;
}

void fmd_min_heap_free(fmd_min_heap_t *heap) {
    free(heap->data);
    free(heap);
}

bool fmd_min_heap_push(fmd_min_heap_t *heap, void *key, void *value) {
    if (heap->size == heap->capacity) {
        heap->data = realloc(heap->data, (heap->capacity *= 2) * sizeof(fmd_min_heap_kv_t));
        if(!heap->data) {
            return false;
        }
    }
    heap->data[heap->size].key = key;
    heap->data[heap->size].value = value;
    fmd_min_heap_sift_up(heap, heap->size);
    heap->size++;
    return true;
}

bool fmd_min_heap_pop(fmd_min_heap_t *heap, void **key, void **value) {
    if (heap->size == 0) {
        *key = NULL;
        *value = NULL;
        return false;
    }
    *key = heap->data[0].key;
    *value = heap->data[0].value;
    heap->data[0] = heap->data[--heap->size];
    fmd_min_heap_sift_down(heap, 0);
    return true;
}

bool fmd_min_heap_peek(fmd_min_heap_t *heap, void **key, void **value) {
    if (heap->size == 0) {
        *key = NULL;
        *value = NULL;
        return false;
    }
    *key = heap->data[0].key;
    *value = heap->data[0].value;
    return true;
}

void fmd_min_heap_sift_up(fmd_min_heap_t *heap, int_t index) {
    while (index > 0) {
        int_t parent = (index-1)>>1;
        if (heap->key_f->comp_f(heap->data[parent].key, heap->data[index].key) <= 0) {
            break;
        }
        fmd_min_heap_kv_t tmp = heap->data[index];
        heap->data[index] = heap->data[parent];
        heap->data[parent] = tmp;
        index = parent;
    }
}

void fmd_min_heap_sift_down(fmd_min_heap_t *heap, int_t index) {
    fmd_min_heap_kv_t tmp = heap->data[index];
    int_t index_2p1 = (index<<1)+1;
    while (index_2p1 < heap->size) {
        int_t smallest = index_2p1;
        int_t right = index_2p1+1;
        if (right < heap->size && heap->key_f->comp_f(heap->data[right].key, heap->data[smallest].key) < 0) {
            smallest = right;
        }
        if (heap->key_f->comp_f(tmp.key, heap->data[smallest].key) <= 0) {
            break;
        }
        heap->data[index] = heap->data[smallest];
        index = smallest;
        index_2p1 = (index<<1)+1;
    }
    heap->data[index] = tmp;
}

void fmd_min_heap_sift_down_old(fmd_min_heap_t *heap, int_t index) {
    while (2 * index + 1 < heap->size) {
        int_t smallest = 2 * index + 1;
        int_t right = 2 * index + 2;
        if (right < heap->size && heap->key_f->comp_f(heap->data[right].key, heap->data[smallest].key) < 0) {
            smallest = right;
        }
        if (heap->key_f->comp_f(heap->data[index].key, heap->data[smallest].key) <= 0) {
            break;
        }
        fmd_min_heap_kv_t tmp = heap->data[index];
        heap->data[index] = heap->data[smallest];
        heap->data[smallest] = tmp;
        index = smallest;
    }
}

int fmd_min_heap_comp(fmd_min_heap_t *h1, fmd_min_heap_t *h2) {
    if (h1->size != h2->size) {
        return (int)(h1->size - h2->size);
    }
    for (int_t i = 0; i < h1->size; ++i) {
        int key_comp = h1->key_f->comp_f(h1->data[i].key, h2->data[i].key);
        if (key_comp != 0) {
            return key_comp;
        }
    }
    return 0;
}

uint_t fmd_min_heap_hash(fmd_min_heap_t *heap) {
    const uint_t prime = 1099511628211LLU;
    uint_t hash = 14695981039346656037LLU;
    for (int_t i = 0; i < heap->size; ++i) {
        uint_t key_hash = heap->key_f->hash_f(heap->data[i].key);
        hash ^= key_hash;
        hash *= prime;
    }
    return hash;
}

fmd_min_heap_t *fmd_min_heap_copy(fmd_min_heap_t *heap) {
    fmd_min_heap_t *copy = (fmd_min_heap_t *)malloc(sizeof(fmd_min_heap_t));
    copy->capacity = heap->capacity;
    copy->size = heap->size;
    copy->key_f = heap->key_f;
    copy->val_f = heap->val_f;
    copy->data = calloc(copy->capacity, sizeof(fmd_min_heap_kv_t));
    for (int_t i = 0; i < copy->size; ++i) {
        copy->data[i] = heap->data[i];
    }
    return copy;
}