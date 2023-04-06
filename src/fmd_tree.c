#include "fmd_tree.h"

fmd_tree_node_t *fmd_tree_node_init(void* key, void* value, fmd_tree_node_t *parent) {
    fmd_tree_node_t *node = (fmd_tree_node_t *)calloc(1, sizeof(fmd_tree_node_t));
    if(!node) {
        return NULL;
    }
    node->key = key;
    node->value = value;
    node->left = NULL;
    node->right = NULL;
    node->parent = parent;
    node->color = RED;
    return node;
}

fmd_tree_node_t *fmd_tree_node_grandparent(fmd_tree_node_t *node) {
    return !node || !node->parent ? NULL : node->parent->parent;
}

void fmd_tree_node_rotate_left(fmd_tree_node_t **root, fmd_tree_node_t *x) {
    fmd_tree_node_t *y = x->right;
    x->right = y->left;
    if (y->left) {
        y->left->parent = x;
    }
    y->parent = x->parent;
    if (!x->parent) {
        (*root) = y;
    } else if (x == x->parent->left) {
        x->parent->left = y;
    } else {
        x->parent->right = y;
    }
    y->left = x;
    x->parent = y;
}

void fmd_tree_node_rotate_right(fmd_tree_node_t **root, fmd_tree_node_t *x) {
    fmd_tree_node_t *y = x->left;
    x->left = y->right;
    if (y->right) {
        y->right->parent = x;
    }
    y->parent = x->parent;
    if (!x->parent) {
        (*root) = y;
    } else if (x == x->parent->right) {
        x->parent->right = y;
    } else {
        x->parent->left = y;
    }
    y->right = x;
    x->parent = y;
}

void fmd_tree_node_fix_violation(fmd_tree_node_t **root, fmd_tree_node_t *node) {
    fmd_tree_node_t *parent = NULL;
    fmd_tree_node_t *grandparent = NULL;
    while ((node != *root) && (node->color == RED) &&
           (node->parent->color == RED)) {
        parent = node->parent;
        grandparent = fmd_tree_node_grandparent(node);
        if (parent == grandparent->left) {
            fmd_tree_node_t *uncle = grandparent->right;
            if (uncle && uncle->color == RED) {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            } else {
                if (node == parent->right) {
                    fmd_tree_node_rotate_left(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                fmd_tree_node_rotate_right(root, grandparent);
                rbt_color_t tempColor = parent->color;
                parent->color = grandparent->color;
                grandparent->color = tempColor;
                node = parent;
            }
        } else {
            fmd_tree_node_t *uncle = grandparent->left;
            if (uncle && uncle->color == RED) {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            } else {
                if (node == parent->left) {
                    fmd_tree_node_rotate_right(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                fmd_tree_node_rotate_left(root, grandparent);
                rbt_color_t tempColor = parent->color;
                parent->color = grandparent->color;
                grandparent->color = tempColor;
                node = parent;
            }
        }
    }
    (*root)->color = BLACK;
}

bool fmd_tree_node_insert(fmd_tree_node_t **root, void* key, void* value, fmd_fstruct_t *key_f) {
    fmd_tree_node_t *node = (*root);
    fmd_tree_node_t *parent = NULL;
    while (node) {
        parent = node;
        int cmpval = key_f->comp_f(key, node->key);
        if (cmpval < 0) {
            node = node->left;
        } else if (cmpval > 0) {
            node = node->right;
        } else {
            return false; // do not allow duplicates
        }
    }
    fmd_tree_node_t *tmp = fmd_tree_node_init(key, value, parent);
    if (!parent) {
        (*root) = tmp;
    } else if (key_f->comp_f(tmp->key, parent->key) < 0) {
        parent->left = tmp;
    } else {
        parent->right = tmp;
    }
    fmd_tree_node_fix_violation(root, tmp);
    return true;
}

fmd_tree_node_t *fmd_tree_node_search(fmd_tree_node_t *root, void *key, fmd_fstruct_t *key_f) {
    while (root) {
        int cmpval = key_f->comp_f(key,root->key);
        if (cmpval < 0) {
            root = root->left;
        } else if (cmpval > 0) {
            root = root->right;
        } else {
            return root;
        }
    }
    return NULL;
}

pos_t fmd_tree_node_height(fmd_tree_node_t *root) {
    return root ? 1 + MAX2(fmd_tree_node_height(root->left), fmd_tree_node_height(root->right)) : -1;
}

void fmd_tree_node_free(fmd_tree_node_t *node, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    if (!node) {
        return;
    }
    if(key_f && key_f->free_f)
        key_f->free_f(node->key);
    if(val_f && val_f->free_f)
        val_f->free_f(node->value);
    fmd_tree_node_free(node->left, key_f, val_f);
    fmd_tree_node_free(node->right, key_f, val_f);
    free(node);
}

fmd_tree_node_t *fmd_tree_node_copy(fmd_tree_node_t *node, fmd_tree_node_t *parent, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    if (!node) {
        return NULL;
    }
    void *key_cpy = key_f->copy_f(node->key);
    void *val_cpy = val_f->copy_f(node->value);
    fmd_tree_node_t *tmp = fmd_tree_node_init(key_cpy, val_cpy, parent);
    if (!tmp) {
        key_f->free_f(key_cpy);
        val_f->free_f(val_cpy);
        return NULL;
    }
    tmp->color = node->color;
    tmp->left = fmd_tree_node_copy(node->left, tmp, key_f, val_f);
    tmp->right = fmd_tree_node_copy(node->right, tmp, key_f, val_f);
    return tmp;
}

int fmd_tree_node_comp(fmd_tree_node_t *t1, fmd_tree_node_t *t2, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    if (t1 == NULL && t2 == NULL) {
        return 0;
    }
    if (t1 == NULL || t2 == NULL) {
        return -1;
    }
    int key_cmp = key_f->comp_f(t1->key, t2->key);
    int val_cmp = val_f->comp_f(t1->value, t2->value);
    if (key_cmp != 0 || val_cmp != 0 || t1->color != t2->color) {
        return -1;
    }
    int left_cmp = fmd_tree_node_comp(t1->left, t2->left, key_f, val_f);
    int right_cmp = fmd_tree_node_comp(t1->right, t2->right, key_f, val_f);
    if (left_cmp || right_cmp ) {
        return -1;
    }
    return 0;
}

upos_t fmd_tree_node_hash(fmd_tree_node_t *node, fmd_fstruct_t *key_f) {
    uint64_t hash = 14695981039346656037LLU;
    if(node) {
        const uint64_t prime = 1099511628211LLU;
        upos_t h = key_f->hash_f(node->key);
        if(node->left) {
            upos_t hl = fmd_tree_node_hash(node->left, key_f);
            h ^= hl;
            h *= prime;
        }
        if(node->right) {
            upos_t hr = fmd_tree_node_hash(node->right, key_f);
            h ^= hr;
            h *= prime;
        }
        return h;
    } else {
        return hash;
    }
}

// traversals
void fmd_tree_node_preorder(fmd_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        f(root->key, root->value, p);
        fmd_tree_node_preorder(root->left,p,f);
        fmd_tree_node_preorder(root->right,p,f);
    }
}

void fmd_tree_node_inorder(fmd_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        fmd_tree_node_inorder(root->left,p,f);
        f(root->key, root->value, p);
        fmd_tree_node_inorder(root->right,p,f);
    }
}

void fmd_tree_node_postorder(fmd_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        fmd_tree_node_postorder(root->left,p,f);
        fmd_tree_node_postorder(root->right,p,f);
        f(root->key, root->value, p);
    }
}

void fmd_tree_init(fmd_tree_t **tree, fmd_fstruct_t *key_f, fmd_fstruct_t *val_f) {
    fmd_tree_t *t = calloc(1, sizeof(fmd_tree_t));
    if(!t) {
        *tree = NULL;
        return;
    }
    t->key_f = key_f;
    t->val_f = val_f;
    t->no_items = 0;
    t->root = NULL;
}

bool fmd_tree_insert(fmd_tree_t *tree, void *key, void *value) {
    bool inserted = fmd_tree_node_insert(&tree->root, key, value, tree->key_f);
    if(inserted) ++tree->no_items;
    return inserted;
}

bool fmd_tree_search(fmd_tree_t *tree, void *key, void **value) {
    fmd_tree_node_t *node = fmd_tree_node_search(tree->root, key, tree->key_f);
    if(!node) return false;
    *value = node->value;
    return true;
}

void fmd_tree_height(fmd_tree_t *tree, pos_t *height) {
    pos_t h = fmd_tree_node_height(tree->root);
    *height = h;
}

void fmd_tree_free(fmd_tree_t *tree) {
    if(tree) {
        fmd_tree_node_free(tree->root, tree->key_f, tree->val_f);
        free(tree);
    }
}

int fmd_tree_comp(fmd_tree_t *t1, fmd_tree_t *t2) {
    if(t1->no_items == t2->no_items) return -1;
    return fmd_tree_node_comp(t1->root,t2->root,t1->key_f,t1->val_f);
}

upos_t fmd_tree_hash(fmd_tree_t *tree) {
    return fmd_tree_node_hash(tree->root, tree->key_f);
}

fmd_tree_t *fmd_tree_copy(fmd_tree_t *tree) {
    fmd_tree_t *copy = calloc(1, sizeof(fmd_tree_t));
    if(!copy) return NULL;
    copy->root = fmd_tree_node_copy(tree->root, NULL, tree->key_f, tree->val_f);
    if(!copy->root) {
        free(copy);
        return NULL;
    }
    copy->key_f = tree->key_f;
    copy->val_f = tree->val_f;
    copy->no_items = tree->no_items;
    return copy;
}

void fmd_tree_preorder(fmd_tree_t *tree, void *p, ftrav_kv f) {
    fmd_tree_node_preorder(tree->root, p, f);
}

void fmd_tree_inorder(fmd_tree_t *tree, void *p, ftrav_kv f) {
    fmd_tree_node_inorder(tree->root, p, f);
}

void fmd_tree_postorder(fmd_tree_t *tree, void *p, ftrav_kv f) {
    fmd_tree_node_postorder(tree->root, p, f);
}