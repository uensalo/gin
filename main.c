#include <stdio.h>
#include <stdlib.h>

typedef enum { RED, BLACK } Color;

typedef struct Node {
    int key;
    struct Node *left;
    struct Node *right;
    struct Node *parent;
    Color color;
} Node;

Node *createNode(int key, Node *parent) {
    Node *node = (Node *)malloc(sizeof(Node));
    node->key = key;
    node->left = NULL;
    node->right = NULL;
    node->parent = parent;
    node->color = RED;
    return node;
}

Node *grandparent(Node *node) {
    if (node == NULL || node->parent == NULL) {
        return NULL;
    }
    return node->parent->parent;
}

Node *uncle(Node *node) {
    Node *g = grandparent(node);
    if (g == NULL) {
        return NULL;
    }
    if (node->parent == g->left) {
        return g->right;
    } else {
        return g->left;
    }
}

void rotateLeft(Node **root, Node *x) {
    Node *y = x->right;
    x->right = y->left;

    if (y->left != NULL) {
        y->left->parent = x;
    }
    y->parent = x->parent;

    if (x->parent == NULL) {
        (*root) = y;
    } else if (x == x->parent->left) {
        x->parent->left = y;
    } else {
        x->parent->right = y;
    }

    y->left = x;
    x->parent = y;
}

void rotateRight(Node **root, Node *x) {
    Node *y = x->left;
    x->left = y->right;

    if (y->right != NULL) {
        y->right->parent = x;
    }
    y->parent = x->parent;

    if (x->parent == NULL) {
        (*root) = y;
    } else if (x == x->parent->right) {
        x->parent->right = y;
    } else {
        x->parent->left = y;
    }

    y->right = x;
    x->parent = y;
}

void fixViolation(Node **root, Node *node) {
    Node *parent = NULL;
    Node *grandParent = NULL;

    while ((node != *root) && (node->color == RED) &&
           (node->parent->color == RED)) {
        parent = node->parent;
        grandParent = grandparent(node);

        if (parent == grandParent->left) {
            Node *uncle = grandParent->right;

            if (uncle != NULL && uncle->color == RED) {
                grandParent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandParent;
            } else {
                if (node == parent->right) {
                    rotateLeft(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                rotateRight(root, grandParent);
                Color tempColor = parent->color;
                parent->color = grandParent->color;
                grandParent->color = tempColor;
                node = parent;
            }
        } else {
            Node *uncle = grandParent->left;

            if (uncle != NULL && uncle->color == RED) {
                grandParent->color = RED;
                parent->color = BLACK;
                uncle->color = BLACK;
                node = grandParent;
            } else {
                if (node == parent->left) {
                    rotateRight(root, parent);
                    node = parent;
                    parent = node->parent;
                }
                rotateLeft(root, grandParent);
                Color tempColor = parent->color;
                parent->color = grandParent->color;
                grandParent->color = tempColor;
                node = parent;
            }
        }
    }

    (*root)->color = BLACK;
}

void insert(Node **root, int key) {
    Node *node = (*root);
    Node *parent = NULL;

    while (node != NULL) {
        parent = node;
        if (key < node->key) {
            node = node->left;
        } else if (key > node->key) {
            node = node->right;
        } else {
            return;
        }
    }

    Node *newNode = createNode(key, parent);
    if (parent == NULL) {
        (*root) = newNode;
    } else if (newNode->key < parent->key) {
        parent->left = newNode;
    } else {
        parent->right = newNode;
    }

    fixViolation(root, newNode);
}

void inorderTraversal(Node *root) {
    if (root == NULL) {
        return;
    }

    inorderTraversal(root->left);
    printf("%d ", root->key);
    inorderTraversal(root->right);
}

Node *search(Node *root, int key) {
    while (root != NULL) {
        if (key < root->key) {
            root = root->left;
        } else if (key > root->key) {
            root = root->right;
        } else {
            return root;
        }
    }
    return NULL;
}

int main() {
    Node *root = NULL;

    int no_elems = 10000;

    for(int i = 1000; i >= 0; i--) {
        insert(&root, i);
    }

    for(int i = 0; i <= 1000; i++) {
        Node* l = search(root, i);
        printf("%d ",l->key);
    }

    printf("Inorder Traversal of the Red-Black tree:\n");
    inorderTraversal(root);
    printf("\n");

    return 0;
}
/*
int main() {
    fmd_bs_t *bs;
    int no_writes = 2022;
    fmd_bs_init_reserve(&bs,1);


    int pointer = 0;
    for (int i = 0; i < no_writes; i++) {
        word_t val = (i % 7) * (i % 11) * (i % 17);
        word_t sz = (i * i) % 13;
        fmd_bs_init_write_word(bs, pointer, val, sz);
        printf("%d ", ~mask_shift_left_64[sz] & val);
        pointer += sz;
    }
    printf("\n");
    pointer = 0;
    for (int i = 0; i < no_writes; i++) {
        word_t val = 0;
        word_t sz = (i * i) % 13;
        fmd_bs_init_read_word(bs, pointer, sz, &val);
        printf("%d ", val);
        pointer += sz;
    }

}
*/