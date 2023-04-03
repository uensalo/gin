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

typedef struct fmd_tree {
    fmd_tree_node_t *root;
    fmd_fstruct_t *key_f;
    fmd_fstruct_t *val_f;
    pos_t no_items;
} fmd_tree_t;

fmd_tree_node_t *fmd_tree_node_init(void* key, void* value, fmd_tree_node_t *parent);
fmd_tree_node_t *fmd_tree_node_grandparent(fmd_tree_node_t *node);
void fmd_tree_node_rotate_left(fmd_tree_node_t **root, fmd_tree_node_t *x);
void fmd_tree_node_rotate_right(fmd_tree_node_t **root, fmd_tree_node_t *x);
void fmd_tree_node_fix_violation(fmd_tree_node_t **root, fmd_tree_node_t *node);
bool fmd_tree_node_insert(fmd_tree_node_t **root, void* key, void* value, fmd_fstruct_t *key_f);
fmd_tree_node_t *fmd_tree_node_search(fmd_tree_node_t *root, void* key, fmd_fstruct_t *key_f);
pos_t fmd_tree_node_height(fmd_tree_node_t *root);
void fmd_tree_node_free(fmd_tree_node_t *node, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);
fmd_tree_node_t *fmd_tree_node_copy(fmd_tree_node_t *node, fmd_tree_node_t *parent, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f);

void fmd_tree_init(fmd_tree_t **tree);
void fmd_tree_insert(fmd_tree_t *tree, void *key, void *value);
void fmd_tree_search(fmd_tree_t *tree, void *key, void **value);
void fmd_tree_height(fmd_tree_t *tree, pos_t *height);
void fmd_tree_free(fmd_tree_t *tree);

#endif //FMD_FMD_TREE_H
