#include "fmd_table.h"

void fmd_table_init(fmd_table_t **table, pos_t capacity, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    fmd_table_t *ht = calloc(1,sizeof(fmd_table_t));
    if(!ht) return;
    ht->size = 0;
    ht->capacity = capacity;
    ht->roots = calloc(capacity,sizeof(fmd_tree_node_t*));
    ht->key_f = key_f;
    ht->val_f = val_f;
    ht->items_per_bucket = calloc(ht->capacity, sizeof(pos_t));
    *table = ht;
}

void fmd_table_free(fmd_table_t *table) {
    if(!table) return;
    for(pos_t i = 0; i < table->capacity; i++) {
        fmd_tree_node_free(table->roots[i], table->key_f, table->val_f);
    }
    free(table->items_per_bucket);
    free(table->roots);
    free(table);
}

bool fmd_table_insert(fmd_table_t *table, void *key, void *value) {
    upos_t index = table->key_f->hash_f(key) % table->capacity;
    bool inserted = fmd_tree_node_insert(&table->roots[index], key, value, table->key_f);
    if(inserted) {
        ++table->size;
        ++table->items_per_bucket[index];
        return true;
    }
    return false;
}

bool fmd_table_lookup(fmd_table_t *table, void *key, void **value) {
    upos_t index = table->key_f->hash_f(key) % table->capacity;
    fmd_tree_node_t *root = table->roots[index];
    fmd_tree_node_t *node = fmd_tree_node_search(root, key, table->key_f);
    if(node) {
        *value = node->value;
        return true;
    }
    *value = NULL;
    return false;
}

void fmd_table_rehash_helper_(fmd_tree_node_t *node, void *nt) {
    fmd_table_t *new_table = (fmd_table_t*)nt;
    fmd_table_insert(new_table, node->key, node->value);
}

void fmd_table_rehash(fmd_table_t **table, pos_t new_capacity) {
    if(!table || !(*table)) return;
    fmd_table_t *ht = *table;
    fmd_table_t *new_ht;
    fmd_table_init(&new_ht, new_capacity, ht->key_f, ht->val_f);
    for (int i = 0; i < ht->capacity; i++) {
        fmd_tree_node_preorder(ht->roots[i], new_ht, fmd_table_rehash_helper_);
    }
    ht->key_f = NULL;
    ht->val_f = NULL;
    fmd_table_free(ht);
    *table = new_ht;
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
    copy->items_per_bucket = calloc(table->capacity, sizeof(pos_t));
    copy->capacity = table->capacity;
    for(pos_t i = 0; i < table->capacity; i++) {
        copy->roots[i] = fmd_tree_node_copy(table->roots[i], NULL, table->key_f,table->val_f);
    }
    copy->size = table->size;
    return copy;
}

upos_t fmd_table_hash(fmd_table_t *table) {
    const upos_t prime = 1099511628211LLU;
    upos_t hash = 14695981039346656037LLU;
    for (pos_t i = 0; i < table->capacity; ++i) {
        upos_t hash_item = fmd_tree_node_hash(table->roots[i], table->key_f);
        hash ^= hash_item;
        hash *= prime;
    }
    return hash;
}