/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_tree.h is part of fmd
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
#ifndef FMD_FMD_TREE_H
#define FMD_FMD_TREE_H
#include "fmd_common.h"

typedef enum { RED, BLACK } rbt_color_t;

typedef struct fmd_tree_node_ {
    void *key;
    void *value;
    struct fmd_tree_node_ *left;
    struct fmd_tree_node_ *right;
    struct fmd_tree_node_ *parent;
    rbt_color_t color;
} fmd_tree_node_t;

fmd_tree_node_t *fmd_tree_node_init(void* key, void* value, fmd_tree_node_t *parent);
fmd_tree_node_t *fmd_tree_node_grandparent(fmd_tree_node_t *node);
void fmd_tree_node_rotate_left(fmd_tree_node_t **root, fmd_tree_node_t *x);
void fmd_tree_node_rotate_right(fmd_tree_node_t **root, fmd_tree_node_t *x);
void fmd_tree_node_fix_violation(fmd_tree_node_t **root, fmd_tree_node_t *node);
bool fmd_tree_node_insert(fmd_tree_node_t **root, void* key, void* value, fmd_fstruct_t *key_f);
fmd_tree_node_t *fmd_tree_node_search(fmd_tree_node_t *root, void* key, fmd_fstruct_t *key_f);
bool fmd_tree_node_replace(fmd_tree_node_t *root, void* key, void* value, void** old_value, fmd_fstruct_t *key_f);
bool fmd_tree_node_replace_if_not_insert(fmd_tree_node_t **root, void* key, void* value, void** old_value, fmd_fstruct_t *key_f);
int_t fmd_tree_node_height(fmd_tree_node_t *root);
void fmd_tree_node_free(fmd_tree_node_t *node, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
fmd_tree_node_t *fmd_tree_node_copy(fmd_tree_node_t *node, fmd_tree_node_t *parent, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
int fmd_tree_node_comp(fmd_tree_node_t *n1, fmd_tree_node_t *n2, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
uint_t fmd_tree_node_hash(fmd_tree_node_t *node, fmd_fstruct_t *key_f);
void fmd_tree_node_preorder(fmd_tree_node_t *root, void *p, ftrav_kv f);
void fmd_tree_node_inorder(fmd_tree_node_t *root, void *p, ftrav_kv f);
void fmd_tree_node_postorder(fmd_tree_node_t *root, void *p, ftrav_kv f);

// clean, exposed API to the mess above
typedef struct fmd_tree {
    fmd_tree_node_t *root;
    fmd_fstruct_t *key_f;
    fmd_fstruct_t *val_f;
    int_t no_items;
} fmd_tree_t;

void fmd_tree_init(fmd_tree_t **tree, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
bool fmd_tree_insert(fmd_tree_t *tree, void *key, void *value);
bool fmd_tree_search(fmd_tree_t *tree, void *key, void **value);
bool fmd_tree_replace(fmd_tree_t *tree, void *key, void *value, void **old_value);
bool fmd_tree_replace_if_not_insert(fmd_tree_t *tree, void *key, void *value, void **old_value);
void fmd_tree_height(fmd_tree_t *tree, int_t *height);
void fmd_tree_free(fmd_tree_t *tree);
int fmd_tree_comp(fmd_tree_t *t1, fmd_tree_t *t2);
uint_t fmd_tree_hash(fmd_tree_t *tree);
fmd_tree_t *fmd_tree_copy(fmd_tree_t *tree);

// traversals
void fmd_tree_preorder(fmd_tree_t *tree, void *p, ftrav_kv f);
void fmd_tree_inorder(fmd_tree_t *tree, void *p, ftrav_kv f);
void fmd_tree_postorder(fmd_tree_t *tree, void *p, ftrav_kv f);

static fmd_fstruct_t fmd_fstruct_tree = {
        (fcomp) fmd_tree_comp,
        (fhash) fmd_tree_hash,
        (ffree) fmd_tree_free,
        (fcopy) fmd_tree_copy
};

#endif //FMD_FMD_TREE_H
