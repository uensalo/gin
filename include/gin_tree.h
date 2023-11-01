/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_tree.h is part of gin
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
#ifndef GIN_GIN_TREE_H
#define GIN_GIN_TREE_H
#include "gin_common.h"

typedef enum { RED, BLACK } rbt_color_t;

typedef struct gin_tree_node_ {
    void *key;
    void *value;
    struct gin_tree_node_ *left;
    struct gin_tree_node_ *right;
    struct gin_tree_node_ *parent;
    rbt_color_t color;
} gin_tree_node_t;

gin_tree_node_t *gin_tree_node_init(void* key, void* value, gin_tree_node_t *parent);
gin_tree_node_t *gin_tree_node_grandparent(gin_tree_node_t *node);
void gin_tree_node_rotate_left(gin_tree_node_t **root, gin_tree_node_t *x);
void gin_tree_node_rotate_right(gin_tree_node_t **root, gin_tree_node_t *x);
void gin_tree_node_fix_violation(gin_tree_node_t **root, gin_tree_node_t *node);
bool gin_tree_node_insert(gin_tree_node_t **root, void* key, void* value, gin_fstruct_t *key_f);
gin_tree_node_t *gin_tree_node_search(gin_tree_node_t *root, void* key, gin_fstruct_t *key_f);
bool gin_tree_node_replace(gin_tree_node_t *root, void* key, void* value, void** old_value, gin_fstruct_t *key_f);
bool gin_tree_node_replace_if_not_insert(gin_tree_node_t **root, void* key, void* value, void** old_value, gin_fstruct_t *key_f);
int_t gin_tree_node_height(gin_tree_node_t *root);
void gin_tree_node_free(gin_tree_node_t *node, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
gin_tree_node_t *gin_tree_node_copy(gin_tree_node_t *node, gin_tree_node_t *parent, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
int gin_tree_node_comp(gin_tree_node_t *n1, gin_tree_node_t *n2, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
uint_t gin_tree_node_hash(gin_tree_node_t *node, gin_fstruct_t *key_f);
void gin_tree_node_preorder(gin_tree_node_t *root, void *p, ftrav_kv f);
void gin_tree_node_inorder(gin_tree_node_t *root, void *p, ftrav_kv f);
void gin_tree_node_postorder(gin_tree_node_t *root, void *p, ftrav_kv f);

// clean, exposed API to the mess above
typedef struct gin_tree {
    gin_tree_node_t *root;
    gin_fstruct_t *key_f;
    gin_fstruct_t *val_f;
    int_t no_items;
} gin_tree_t;

void gin_tree_init(gin_tree_t **tree, gin_fstruct_t *key_f, gin_fstruct_t *val_f);
bool gin_tree_insert(gin_tree_t *tree, void *key, void *value);
bool gin_tree_search(gin_tree_t *tree, void *key, void **value);
bool gin_tree_replace(gin_tree_t *tree, void *key, void *value, void **old_value);
bool gin_tree_replace_if_not_insert(gin_tree_t *tree, void *key, void *value, void **old_value);
void gin_tree_height(gin_tree_t *tree, int_t *height);
void gin_tree_free(gin_tree_t *tree);
int gin_tree_comp(gin_tree_t *t1, gin_tree_t *t2);
uint_t gin_tree_hash(gin_tree_t *tree);
gin_tree_t *gin_tree_copy(gin_tree_t *tree);

// traversals
void gin_tree_preorder(gin_tree_t *tree, void *p, ftrav_kv f);
void gin_tree_inorder(gin_tree_t *tree, void *p, ftrav_kv f);
void gin_tree_postorder(gin_tree_t *tree, void *p, ftrav_kv f);

static gin_fstruct_t gin_fstruct_tree = {
        (fcomp) gin_tree_comp,
        (fhash) gin_tree_hash,
        (ffree) gin_tree_free,
        (fcopy) gin_tree_copy
};

#endif //GIN_GIN_TREE_H
