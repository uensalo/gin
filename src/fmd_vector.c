#include "fmd_vector.h"

void fmd_vector_init(fmd_vector_t **vec, pos_t initial_capacity) {
    fmd_vector_t *v = calloc(1, sizeof(fmd_vector_t));
    if(!v) {
        *vec = NULL;
        return;
    }
    v->size = 0;
    v->capacity = initial_capacity;
    v->data = malloc(initial_capacity * sizeof(void *));
    *vec = v;
}

void fmd_vector_free(fmd_vector_t *vec) {
    free(vec->data);
    free(vec);
}

void fmd_vector_grow(fmd_vector_t *vec) {
    vec->capacity *= 2;
    vec->data = realloc(vec->data, vec->capacity * sizeof(void *));
}

void fmd_vector_shrink(fmd_vector_t *vec) {
    vec->capacity /= 2;
    vec->data = realloc(vec->data, vec->capacity * sizeof(void *));
}

void fmd_vector_append(fmd_vector_t *vec, void *value) {
    if (vec->size == vec->capacity) {
        fmd_vector_grow(vec);
    }
    vec->data[vec->size++] = value;
}

void fmd_vector_pop(fmd_vector_t *vec, void **item) {
    void *tmp = vec->data[vec->size-1];
    vec->size--;
    if (vec->size < vec->capacity / 4) {
        fmd_vector_shrink(vec);
    }
    if(item)
        *item = tmp;
}

void fmd_vector_insert(fmd_vector_t *vec, pos_t index, void *value) {
    if (vec->size == vec->capacity) {
        fmd_vector_grow(vec);
    }
    memmove(&vec->data[index + 1], &vec->data[index], (vec->size - index) * sizeof(void *));
    vec->data[index] = value;
    vec->size++;
}

void fmd_vector_delete(fmd_vector_t *vec, pos_t index, void **item) {
    void *tmp = vec->data[index];
    memmove(&vec->data[index], &vec->data[index + 1], (vec->size - index - 1) * sizeof(void *));
    vec->size--;
    if (vec->size < vec->capacity / 4) {
        fmd_vector_shrink(vec);
    }
    if(item)
        *item = tmp;
}

void fmd_vector_swap_(void **a, void **b) {
    void *temp = *a;
    *a = *b;
    *b = temp;
}

void fmd_qs_helper_(void **arr, pos_t left, pos_t right, fcomp comp_f) {
    if (left < right) {
        if (right - left <= 10) {
            for (pos_t i = left + 1; i <= right; i++) {
                void *key = arr[i];
                pos_t j = i - 1;

                while (j >= left && comp_f(arr[j], key) > 0) {
                    arr[j + 1] = arr[j];
                    j--;
                }
                arr[j + 1] = key;
            }
        } else {
            pos_t pivot_index = left + rand() % (right - left + 1);
            void *pivot = arr[pivot_index];
            fmd_vector_swap_(&arr[pivot_index], &arr[left]);
            pos_t i = left + 1;
            pos_t j = right;
            while (true) {
                while (i <= j && comp_f(arr[i], pivot) < 0) {
                    i++;
                }
                while (i <= j && comp_f(arr[j], pivot) > 0) {
                    j--;
                }
                if (i <= j) {
                    fmd_vector_swap_(&arr[i], &arr[j]);
                    i++;
                    j--;
                } else {
                    break;
                }
            }
            fmd_vector_swap_(&arr[left], &arr[j]);
            pivot_index = j;
            fmd_qs_helper_(arr, left, pivot_index, comp_f);
            fmd_qs_helper_(arr, pivot_index + 1, right, comp_f);
        }
    }
}

void fmd_vector_sort(fmd_vector_t *vec, fcomp comp_f) {
    fmd_qs_helper_(vec->data, 0, vec->size - 1, comp_f);
}