#include "fmd_table.h"

void fmd_table_init(fmd_table_t **table, int_t capacity, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    fmd_table_t *ht = calloc(1,sizeof(fmd_table_t));
    if(!ht) return;
    ht->size = 0;
    ht->capacity = capacity;
    ht->roots = calloc(capacity,sizeof(fmd_tree_node_t*));
    ht->key_f = key_f;
    ht->val_f = val_f;
    ht->items_per_bucket = calloc(ht->capacity, sizeof(int_t));
    *table = ht;
}

void fmd_table_free(fmd_table_t *table) {
    if(!table) return;
    for(int_t i = 0; i < table->capacity; i++) {
        fmd_tree_node_free(table->roots[i], table->key_f, table->val_f);
    }
    free(table->items_per_bucket);
    free(table->roots);
    free(table);
}

bool fmd_table_insert(fmd_table_t *table, void *key, void *value) {
    uint_t index = table->key_f->hash_f(key) % table->capacity;
    bool inserted = fmd_tree_node_insert(&table->roots[index], key, value, table->key_f);
    if(inserted) {
        ++table->size;
        ++table->items_per_bucket[index];
        if((double)table->items_per_bucket[index] > (double)table->capacity * FMD_REHASH_FACTOR) {
            int_t next_size = fmd_next_prime_after(2 * table->capacity);
            fmd_table_rehash(&table, next_size);
        }
        return true;
    }
    return false;
}

bool fmd_table_lookup(fmd_table_t *table, void *key, void **value) {
    uint_t index = table->key_f->hash_f(key) % table->capacity;
    fmd_tree_node_t *root = table->roots[index];
    fmd_tree_node_t *node = fmd_tree_node_search(root, key, table->key_f);
    if(node) {
        *value = node->value;
        return true;
    }
    *value = NULL;
    return false;
}

void fmd_table_rehash_helper_(void* key, void* value, void *nt) {
    typedef struct rehash_params_ {
        fmd_tree_node_t **buckets;
        int_t *items_per_bucket;
        fmd_fstruct_t *key_f;
        fmd_fstruct_t *val_f;
        int_t new_cap;
    } rehash_params_t;

    rehash_params_t *p = (rehash_params_t*)nt;
    uint_t index = p->key_f->hash_f(key) % p->new_cap;
    bool inserted = fmd_tree_node_insert(&p->buckets[index], key, value, p->key_f);
    ++p->items_per_bucket[index];
}

void fmd_table_rehash(fmd_table_t **table, int_t new_capacity) {
    if(!table || !(*table)) return;
    fmd_table_t *ht = *table;
    fmd_tree_node_t **buckets = calloc(new_capacity, sizeof(fmd_tree_node_t*));
    int_t *items_per_bucket = calloc(new_capacity, sizeof(int_t));
    typedef struct rehash_params_ {
        fmd_tree_node_t **buckets;
        int_t *items_per_bucket;
        fmd_fstruct_t *key_f;
        fmd_fstruct_t *val_f;
        int_t new_cap;
    } rehash_params_t;
    rehash_params_t p;
    p.buckets = buckets;
    p.items_per_bucket = items_per_bucket;
    p.key_f = ht->key_f;
    p.val_f = ht->val_f;
    p.new_cap = new_capacity;

    for (int i = 0; i < ht->capacity; i++) {
        fmd_tree_node_preorder(ht->roots[i], &p, (ftrav_kv) fmd_table_rehash_helper_);
        fmd_tree_node_free(ht->roots[i], NULL, NULL);
    }
    free(ht->roots);
    free(ht->items_per_bucket);

    ht->items_per_bucket = items_per_bucket;
    ht->roots = buckets;
    ht->capacity = new_capacity;
}

// todo: traverse one table and call lookup on the other
// boring to write, and not needed unless nested tables are a thing
int fmd_table_comp(fmd_table_t *t1, fmd_table_t *t2) {
    return -1;
}

fmd_table_t *fmd_table_copy(fmd_table_t *table) {
    fmd_table_t *copy = calloc(1, sizeof(fmd_table_t));
    copy->key_f = table->key_f;
    copy->val_f = table->val_f;
    copy->roots = calloc(table->capacity, sizeof(fmd_tree_node_t*));
    copy->items_per_bucket = calloc(table->capacity, sizeof(int_t));
    copy->capacity = table->capacity;
    for(int_t i = 0; i < table->capacity; i++) {
        copy->roots[i] = fmd_tree_node_copy(table->roots[i], NULL, table->key_f,table->val_f);
    }
    copy->size = table->size;
    return copy;
}

uint_t fmd_table_hash(fmd_table_t *table) {
    const uint_t prime = 1099511628211LLU;
    uint_t hash = 14695981039346656037LLU;
    for (int_t i = 0; i < table->capacity; ++i) {
        uint_t hash_item = fmd_tree_node_hash(table->roots[i], table->key_f);
        hash ^= hash_item;
        hash *= prime;
    }
    return hash;
}

void fmd_table_traverse(fmd_table_t *table, void *p, ftrav_kv f) {
    if (!table) return;
    for (int_t i = 0; i < table->capacity; i++) {
        fmd_tree_node_preorder(table->roots[i], p, f);
    }
}