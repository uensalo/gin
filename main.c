#include <stdio.h>
#include <stdlib.h>
#include <fmd_tree.h>
#include <time.h>

int main() {
    fmd_tree_node_t *root = NULL;

    int no_elems = 100;

    for(int i = no_elems; i >= 0; i--) {
        fmd_tree_node_insert(&root, i, i*i, &prm_fstruct);
    }

    for(int i = 0; i <= no_elems; i++) {
        //fmd_tree_node_t* l = fmd_tree_node_search(root, i, &prm_fstruct);
        //printf("%d, %d\n",l->key, l->value);
    }
    //printf("%d\n", fmd_tree_node_height(root));
    //fmd_tree_node_free(root, &prm_fstruct, &prm_fstruct);

    fmd_tree_node_t *cpy = fmd_tree_node_copy(root, NULL, &prm_fstruct, &prm_fstruct);

    for(int i = 0; i <= no_elems; i++) {
        fmd_tree_node_t* l = fmd_tree_node_search(cpy, i, &prm_fstruct);
        printf("%d, %d\n",l->key, l->value);
    }

    return 0;
}