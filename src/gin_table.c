/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_table.c is part of gin
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
#include "gin_table.h"

void gin_table_init(gin_table_t **table, int_t capacity, gin_fstruct_t *key_f, gin_fstruct_t *val_f) {
    gin_table_t *ht = calloc(1,sizeof(gin_table_t));
    if(!ht) return;
    ht->size = 0;
    ht->capacity = capacity;
    ht->roots = calloc(capacity,sizeof(gin_tree_node_t*));
    ht->key_f = key_f;
    ht->val_f = val_f;
    ht->items_per_bucket = calloc(ht->capacity, sizeof(int_t));
    *table = ht;
}

void gin_table_free(gin_table_t *table) {
    if(!table) return;
    for(int_t i = 0; i < table->capacity; i++) {
        gin_tree_node_free(table->roots[i], table->key_f, table->val_f);
    }
    free(table->items_per_bucket);
    free(table->roots);
    free(table);
}

bool gin_table_insert(gin_table_t *table, void *key, void *value) {
    uint_t index = table->key_f->hash_f(key) % table->capacity;
    bool inserted = gin_tree_node_insert(&table->roots[index], key, value, table->key_f);
    if(inserted) {
        ++table->size;
        ++table->items_per_bucket[index];
        if((double)table->items_per_bucket[index] > (double)table->capacity * GIN_REHASH_FACTOR) {
            int_t next_size = gin_next_prime_after(2 * table->capacity);
            gin_table_rehash(&table, next_size);
        }
        return true;
    }
    return false;
}

bool gin_table_lookup(gin_table_t *table, void *key, void **value) {
    uint_t index = table->key_f->hash_f(key) % table->capacity;
    gin_tree_node_t *root = table->roots[index];
    gin_tree_node_t *node = gin_tree_node_search(root, key, table->key_f);
    if(node) {
        *value = node->value;
        return true;
    }
    *value = NULL;
    return false;
}

void gin_table_rehash_helper_(void* key, void* value, void *nt) {
    typedef struct rehash_params_ {
        gin_tree_node_t **buckets;
        int_t *items_per_bucket;
        gin_fstruct_t *key_f;
        gin_fstruct_t *val_f;
        int_t new_cap;
    } rehash_params_t;

    rehash_params_t *p = (rehash_params_t*)nt;
    uint_t index = p->key_f->hash_f(key) % p->new_cap;
    bool inserted = gin_tree_node_insert(&p->buckets[index], key, value, p->key_f);
    ++p->items_per_bucket[index];
}

void gin_table_rehash(gin_table_t **table, int_t new_capacity) {
    if(!table || !(*table)) return;
    gin_table_t *ht = *table;
    gin_tree_node_t **buckets = calloc(new_capacity, sizeof(gin_tree_node_t*));
    int_t *items_per_bucket = calloc(new_capacity, sizeof(int_t));
    typedef struct rehash_params_ {
        gin_tree_node_t **buckets;
        int_t *items_per_bucket;
        gin_fstruct_t *key_f;
        gin_fstruct_t *val_f;
        int_t new_cap;
    } rehash_params_t;
    rehash_params_t p;
    p.buckets = buckets;
    p.items_per_bucket = items_per_bucket;
    p.key_f = ht->key_f;
    p.val_f = ht->val_f;
    p.new_cap = new_capacity;

    for (int i = 0; i < ht->capacity; i++) {
        gin_tree_node_preorder(ht->roots[i], &p, (ftrav_kv) gin_table_rehash_helper_);
        gin_tree_node_free(ht->roots[i], NULL, NULL);
    }
    free(ht->roots);
    free(ht->items_per_bucket);

    ht->items_per_bucket = items_per_bucket;
    ht->roots = buckets;
    ht->capacity = new_capacity;
}

// todo: traverse one table and call lookup on the other
// boring to write, and not needed unless nested tables are a thing
int gin_table_comp(gin_table_t *t1, gin_table_t *t2) {
    return -1;
}

gin_table_t *gin_table_copy(gin_table_t *table) {
    gin_table_t *copy = calloc(1, sizeof(gin_table_t));
    copy->key_f = table->key_f;
    copy->val_f = table->val_f;
    copy->roots = calloc(table->capacity, sizeof(gin_tree_node_t*));
    copy->items_per_bucket = calloc(table->capacity, sizeof(int_t));
    copy->capacity = table->capacity;
    for(int_t i = 0; i < table->capacity; i++) {
        copy->roots[i] = gin_tree_node_copy(table->roots[i], NULL, table->key_f,table->val_f);
    }
    copy->size = table->size;
    return copy;
}

uint_t gin_table_hash(gin_table_t *table) {
    const uint_t prime = 1099511628211LLU;
    uint_t hash = 14695981039346656037LLU;
    for (int_t i = 0; i < table->capacity; ++i) {
        uint_t hash_item = gin_tree_node_hash(table->roots[i], table->key_f);
        hash ^= hash_item;
        hash *= prime;
    }
    return hash;
}

void gin_table_traverse(gin_table_t *table, void *p, ftrav_kv f) {
    if (!table) return;
    for (int_t i = 0; i < table->capacity; i++) {
        gin_tree_node_preorder(table->roots[i], p, f);
    }
}