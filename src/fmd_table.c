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
    pos_t index = table->key_f->hash_f(key) % table->capacity;
    fmd_tree_node_t *root = table->roots[index];
    bool inserted = fmd_tree_node_insert(&root, key, value, table->key_f);
    if(inserted) {
        table->size++;
        table->items_per_bucket[index];
        return true;
    }
    return false;
}

void fmd_table_lookup(fmd_table_t *table, void *key, void **value) {
    pos_t index = table->key_f->hash_f(key) % table->capacity;
    fmd_tree_node_t *root = table->roots[index];
    fmd_tree_node_t *node = fmd_tree_node_search(root, key, table->key_f);
    *value = node ? node->value : NULL;
}

void fmd_table_rehash(fmd_table_t **table, pos_t new_capacity) {
    if(!table || !(*table)) return;
    fmd_table_t *ht = *table;
    fmd_table_t *new_ht;
    fmd_table_init(&new_ht, new_capacity, ht->key_f, ht->val_f);

    for (int i = 0; i < ht->capacity; i++) {
        // todo: pre-order-traverse the tree and insert keys one by one to the new one. No copies.
    }

    fmd_table_free(ht);
    *ht = *new_ht;
    free(new_ht);
}