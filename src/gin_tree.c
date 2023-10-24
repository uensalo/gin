/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_tree.c is part of gin
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
#include "gin_tree.h"

gin_tree_node_t *gin_tree_node_init(void* key, void* value, gin_tree_node_t *parent) {
    gin_tree_node_t *node = (gin_tree_node_t *)calloc(1, sizeof(gin_tree_node_t));
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

gin_tree_node_t *gin_tree_node_grandparent(gin_tree_node_t *node) {
    return !node || !node->parent ? NULL : node->parent->parent;
}

void gin_tree_node_rotate_left(gin_tree_node_t **root, gin_tree_node_t *x) {
    gin_tree_node_t *y = x->right;
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

void gin_tree_node_rotate_right(gin_tree_node_t **root, gin_tree_node_t *x) {
    gin_tree_node_t *y = x->left;
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

void gin_tree_node_fix_violation(gin_tree_node_t **root, gin_tree_node_t *node) {
    gin_tree_node_t *parent = NULL;
    gin_tree_node_t *grandparent = NULL;
    while ((node != *root) && (node->color == RED) &&
           (node->parent->color == RED)) {
        parent = node->parent;
        grandparent = gin_tree_node_grandparent(node);
        if (parent == grandparent->left) {
            gin_tree_node_t *uncle = grandparent->right;
            if (uncle && uncle->color == RED) {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            } else {
                if (node == parent->right) {
                    gin_tree_node_rotate_left(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                gin_tree_node_rotate_right(root, grandparent);
                rbt_color_t tempColor = parent->color;
                parent->color = grandparent->color;
                grandparent->color = tempColor;
                node = parent;
            }
        } else {
            gin_tree_node_t *uncle = grandparent->left;
            if (uncle && uncle->color == RED) {
                grandparent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandparent;
            } else {
                if (node == parent->left) {
                    gin_tree_node_rotate_right(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                gin_tree_node_rotate_left(root, grandparent);
                rbt_color_t tempColor = parent->color;
                parent->color = grandparent->color;
                grandparent->color = tempColor;
                node = parent;
            }
        }
    }
    (*root)->color = BLACK;
}

bool gin_tree_node_insert(gin_tree_node_t **root, void* key, void* value, gin_fstruct_t *key_f) {
    gin_tree_node_t *node = (*root);
    gin_tree_node_t *parent = NULL;
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
    gin_tree_node_t *tmp = gin_tree_node_init(key, value, parent);
    if (!parent) {
        (*root) = tmp;
    } else if (key_f->comp_f(tmp->key, parent->key) < 0) {
        parent->left = tmp;
    } else {
        parent->right = tmp;
    }
    gin_tree_node_fix_violation(root, tmp);
    return true;
}

gin_tree_node_t *gin_tree_node_search(gin_tree_node_t *root, void *key, gin_fstruct_t *key_f) {
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

bool gin_tree_node_replace(gin_tree_node_t *root, void* key, void* value, void **old_value, gin_fstruct_t *key_f) {
    while (root) {
        int cmpval = key_f->comp_f(key,root->key);
        if (cmpval < 0) {
            root = root->left;
        } else if (cmpval > 0) {
            root = root->right;
        } else {
            *old_value = root->value;
            root->value = value;
            return true;
        }
    }
    return false;
}

bool gin_tree_node_replace_if_not_insert(gin_tree_node_t **root, void* key, void* value, void** old_value, gin_fstruct_t *key_f) {
    gin_tree_node_t *node = (*root);
    gin_tree_node_t *parent = NULL;
    while (node) {
        parent = node;
        int cmpval = key_f->comp_f(key, node->key);
        if (cmpval < 0) {
            node = node->left;
        } else if (cmpval > 0) {
            node = node->right;
        } else {
            *old_value = node->value;
            node->value = value;
            return true;
        }
    }
    gin_tree_node_t *tmp = gin_tree_node_init(key, value, parent);
    if (!parent) {
        (*root) = tmp;
    } else if (key_f->comp_f(tmp->key, parent->key) < 0) {
        parent->left = tmp;
    } else {
        parent->right = tmp;
    }
    *old_value = NULL;
    gin_tree_node_fix_violation(root, tmp);
    return false;
}

int_t gin_tree_node_height(gin_tree_node_t *root) {
    return root ? 1 + MAX2(gin_tree_node_height(root->left), gin_tree_node_height(root->right)) : -1;
}

void gin_tree_node_free(gin_tree_node_t *node, gin_fstruct_t *key_f, gin_fstruct_t *val_f) {
    if (!node) {
        return;
    }
    if(key_f && key_f->free_f)
        key_f->free_f(node->key);
    if(val_f && val_f->free_f)
        val_f->free_f(node->value);
    gin_tree_node_free(node->left, key_f, val_f);
    gin_tree_node_free(node->right, key_f, val_f);
    free(node);
}

gin_tree_node_t *gin_tree_node_copy(gin_tree_node_t *node, gin_tree_node_t *parent, gin_fstruct_t *key_f, gin_fstruct_t *val_f) {
    if (!node) {
        return NULL;
    }
    void *key_cpy = key_f->copy_f(node->key);
    void *val_cpy = val_f->copy_f(node->value);
    gin_tree_node_t *tmp = gin_tree_node_init(key_cpy, val_cpy, parent);
    if (!tmp) {
        key_f->free_f(key_cpy);
        val_f->free_f(val_cpy);
        return NULL;
    }
    tmp->color = node->color;
    tmp->left = gin_tree_node_copy(node->left, tmp, key_f, val_f);
    tmp->right = gin_tree_node_copy(node->right, tmp, key_f, val_f);
    return tmp;
}

int gin_tree_node_comp(gin_tree_node_t *t1, gin_tree_node_t *t2, gin_fstruct_t *key_f, gin_fstruct_t *val_f) {
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
    int left_cmp = gin_tree_node_comp(t1->left, t2->left, key_f, val_f);
    int right_cmp = gin_tree_node_comp(t1->right, t2->right, key_f, val_f);
    if (left_cmp || right_cmp ) {
        return -1;
    }
    return 0;
}

uint_t gin_tree_node_hash(gin_tree_node_t *node, gin_fstruct_t *key_f) {
    uint64_t hash = 14695981039346656037LLU;
    if(node) {
        const uint64_t prime = 1099511628211LLU;
        uint_t h = key_f->hash_f(node->key);
        if(node->left) {
            uint_t hl = gin_tree_node_hash(node->left, key_f);
            h ^= hl;
            h *= prime;
        }
        if(node->right) {
            uint_t hr = gin_tree_node_hash(node->right, key_f);
            h ^= hr;
            h *= prime;
        }
        return h;
    } else {
        return hash;
    }
}

// traversals
void gin_tree_node_preorder(gin_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        f(root->key, root->value, p);
        gin_tree_node_preorder(root->left,p,f);
        gin_tree_node_preorder(root->right,p,f);
    }
}

void gin_tree_node_inorder(gin_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        gin_tree_node_inorder(root->left,p,f);
        f(root->key, root->value, p);
        gin_tree_node_inorder(root->right,p,f);
    }
}

void gin_tree_node_postorder(gin_tree_node_t *root, void *p, ftrav_kv f) {
    if(root) {
        gin_tree_node_postorder(root->left,p,f);
        gin_tree_node_postorder(root->right,p,f);
        f(root->key, root->value, p);
    }
}

void gin_tree_init(gin_tree_t **tree, gin_fstruct_t *key_f, gin_fstruct_t *val_f) {
    gin_tree_t *t = calloc(1, sizeof(gin_tree_t));
    if(!t) {
        *tree = NULL;
        return;
    }
    t->key_f = key_f;
    t->val_f = val_f;
    t->no_items = 0;
    t->root = NULL;
    *tree = t;
}

bool gin_tree_insert(gin_tree_t *tree, void *key, void *value) {
    bool inserted = gin_tree_node_insert(&tree->root, key, value, tree->key_f);
    if(inserted) ++tree->no_items;
    return inserted;
}

bool gin_tree_search(gin_tree_t *tree, void *key, void **value) {
    gin_tree_node_t *node = gin_tree_node_search(tree->root, key, tree->key_f);
    if(!node) return false;
    *value = node->value;
    return true;
}

bool gin_tree_replace(gin_tree_t *tree, void *key, void *value, void **old_value) {
    void *tmp = NULL;
    bool replaced = gin_tree_node_replace(tree->root, key, value, &tmp, tree->key_f);
    *old_value = tmp;
    return replaced;
}


bool gin_tree_replace_if_not_insert(gin_tree_t *tree, void *key, void *value, void **old_value) {
    void *tmp = NULL;
    bool replaced = gin_tree_node_replace_if_not_insert(&tree->root, key, value, &tmp, tree->key_f);
    if(!replaced) ++tree->no_items;
    *old_value = tmp;
    return replaced;
}

void gin_tree_height(gin_tree_t *tree, int_t *height) {
    int_t h = gin_tree_node_height(tree->root);
    *height = h;
}

void gin_tree_free(gin_tree_t *tree) {
    if(tree) {
        gin_tree_node_free(tree->root, tree->key_f, tree->val_f);
        free(tree);
    }
}

int gin_tree_comp(gin_tree_t *t1, gin_tree_t *t2) {
    if(t1->no_items == t2->no_items) return -1;
    return gin_tree_node_comp(t1->root,t2->root,t1->key_f,t1->val_f);
}

uint_t gin_tree_hash(gin_tree_t *tree) {
    return gin_tree_node_hash(tree->root, tree->key_f);
}

gin_tree_t *gin_tree_copy(gin_tree_t *tree) {
    gin_tree_t *copy = calloc(1, sizeof(gin_tree_t));
    if(!copy) return NULL;
    copy->root = gin_tree_node_copy(tree->root, NULL, tree->key_f, tree->val_f);
    if(!copy->root) {
        free(copy);
        return NULL;
    }
    copy->key_f = tree->key_f;
    copy->val_f = tree->val_f;
    copy->no_items = tree->no_items;
    return copy;
}

void gin_tree_preorder(gin_tree_t *tree, void *p, ftrav_kv f) {
    gin_tree_node_preorder(tree->root, p, f);
}

void gin_tree_inorder(gin_tree_t *tree, void *p, ftrav_kv f) {
    gin_tree_node_inorder(tree->root, p, f);
}

void gin_tree_postorder(gin_tree_t *tree, void *p, ftrav_kv f) {
    gin_tree_node_postorder(tree->root, p, f);
}