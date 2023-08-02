#include "fmd_vector.h"

void fmd_vector_init(fmd_vector_t **vec, int_t initial_capacity, fmd_fstruct_t *f) {
    fmd_vector_t *v = calloc(1, sizeof(fmd_vector_t));
    if(!v) {
        *vec = NULL;
        return;
    }
    v->size = 0;
    v->capacity = initial_capacity;
    v->data = calloc(initial_capacity,sizeof(void *));
    v->f = f;
    *vec = v;
}

void fmd_vector_free(fmd_vector_t *vec) {
    if(vec) {
        if(vec->f && vec->f->free_f && vec->f->free_f != prm_free) {
            for (int_t i = 0; i < vec->size; i++) {
                vec->f->free_f(vec->data[i]);
            }
        }
        free(vec->data);
        free(vec);
    }
}

void fmd_vector_free_disown(fmd_vector_t *vec) {
    if(vec) {
        free(vec->data);
        free(vec);
    }
}

void fmd_vector_grow(fmd_vector_t *vec) {
    vec->data = realloc(vec->data, (vec->capacity *= 2) * sizeof(void *));
}

void fmd_vector_shrink(fmd_vector_t *vec) {
    vec->data = realloc(vec->data, (vec->capacity /= 2) * sizeof(void *));
}

void fmd_vector_fit(fmd_vector_t *vec) {
    vec->data = realloc(vec->data, (vec->capacity = vec->size) * sizeof(void *));
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

void fmd_vector_insert(fmd_vector_t *vec, int_t index, void *value) {
    if (vec->size == vec->capacity) {
        fmd_vector_grow(vec);
    }
    memmove(&vec->data[index + 1], &vec->data[index], (vec->size - index) * sizeof(void*));
    vec->data[index] = value;
    vec->size++;
}

void fmd_vector_delete(fmd_vector_t *vec, int_t index, void **item) {
    void *tmp = vec->data[index];
    memmove(&vec->data[index], &vec->data[index + 1], (vec->size - index - 1) * sizeof(void*));
    vec->size--;
    if (vec->size < vec->capacity / 4) {
        fmd_vector_shrink(vec);
    }
    if(item)
        *item = tmp;
}

void fmd_qs_helper_(void **arr, int_t left, int_t right, fcomp comp_f) {
    if (left < right) {
        if (right - left <= 10) {
            for (int_t i = left + 1; i <= right; i++) {
                void *key = arr[i];
                int_t j = i - 1;

                while (j >= left && comp_f(arr[j], key) > 0) {
                    arr[j + 1] = arr[j];
                    j--;
                }
                arr[j + 1] = key;
            }
        } else {
            int_t pivot_index = left + (int_t)rand() % (right - left + 1);
            void *pivot = arr[pivot_index];
            fmd_swap(&arr[pivot_index], &arr[left]);
            int_t i = left + 1;
            int_t j = right;
            while (true) {
                while (i <= j && comp_f(arr[i], pivot) < 0) {
                    i++;
                }
                while (i <= j && comp_f(arr[j], pivot) > 0) {
                    j--;
                }
                if (i <= j) {
                    fmd_swap(&arr[i], &arr[j]);
                    i++;j--;
                } else {
                    break;
                }
            }
            fmd_swap(&arr[left], &arr[j]);
            pivot_index = j;
            fmd_qs_helper_(arr, left, pivot_index, comp_f);
            fmd_qs_helper_(arr, pivot_index + 1, right, comp_f);
        }
    }
}

void fmd_vector_sort(fmd_vector_t *vec) {
    fmd_qs_helper_(vec->data, 0, vec->size - 1, vec->f->comp_f);
}

void fmd_qs_arg_helper_(void **arr, void **args, int_t left, int_t right, fcomp comp_f) {
    if (left < right) {
        if (right - left <= 10) {
            for (int_t i = left + 1; i <= right; i++) {
                void *key = arr[i];
                void *argkey = args[i];
                int_t j = i - 1;

                while (j >= left && comp_f(arr[j], key) > 0) {
                    arr[j + 1] = arr[j];
                    args[j + 1] = args[j];
                    j--;
                }
                arr[j + 1] = key;
                args[j + 1] = argkey;
            }
        } else {
            int_t pivot_index = left + (int_t)rand() % (right - left + 1);
            void *pivot = arr[pivot_index];
            fmd_swap(&arr[pivot_index], &arr[left]);
            fmd_swap(&args[pivot_index], &args[left]);
            int_t i = left + 1;
            int_t j = right;
            while (true) {
                while (i <= j && comp_f(arr[i], pivot) < 0) {
                    i++;
                }
                while (i <= j && comp_f(arr[j], pivot) > 0) {
                    j--;
                }
                if (i <= j) {
                    fmd_swap(&arr[i], &arr[j]);
                    fmd_swap(&args[i], &args[j]);
                    i++;j--;
                } else {
                    break;
                }
            }
            fmd_swap(&arr[left], &arr[j]);
            fmd_swap(&args[left], &args[j]);
            pivot_index = j;
            fmd_qs_arg_helper_(arr, args, left, pivot_index, comp_f);
            fmd_qs_arg_helper_(arr, args, pivot_index + 1, right, comp_f);
        }
    }
}

void fmd_vector_argsort(fmd_vector_t **args_sorted, fmd_vector_t *vec) {
    fmd_vector_t *args;
    fmd_vector_init(&args, vec->size, &prm_fstruct);
    if(!args) {
        *args_sorted = NULL;
        return;
    }
    for(int_t i = 0; i < vec->size; i++) {
        fmd_vector_append(args, (void*)i);
    }
    fmd_qs_arg_helper_(vec->data, args->data, 0, vec->size - 1, vec->f->comp_f);
    *args_sorted = args;
}

int fmd_vector_comp(fmd_vector_t *v1, fmd_vector_t *v2) {
    if(v1 == v2) return 0;
    if(!v1 && v2 || v1 && !v2) return -1;
    if(v1->size!=v2->size) return -1;
    for(int_t i = 0; i < v1->size; i++) {
        int cmp = v1->f->comp_f(v1->data[i],v2->data[i]);
        if(cmp) return -1;
    }
    return 0;
}

uint_t fmd_vector_hash(fmd_vector_t *v) {
    // fnv again
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;
    for (size_t i = 0; i < v->size; ++i) {
        uint_t hash_item = v->f->hash_f(v->data[i]);
        hash ^= hash_item;
        hash *= prime;
    }
    return hash;
}

fmd_vector_t *fmd_vector_copy(fmd_vector_t *v) {
    fmd_vector_t *cpy = calloc(1, sizeof(fmd_vector_t));
    cpy->data = calloc(v->capacity, sizeof(void*));
    cpy->size = v->size;
    cpy->capacity = v->capacity;
    cpy->f = v->f;
    for(int_t i = 0; i < v->size; i++) {
        cpy->data[i] = v->f->copy_f(v->data[i]);
    }
    return cpy;
}