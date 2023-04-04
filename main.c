#include <stdio.h>
#include <stdlib.h>
#include <fmd_table.h>
#include <fmd_bitstream.h>
#include <time.h>

void print_tree_trav(fmd_tree_node_t *node, void *p) {
    printf("%d, %d\n", node->key, node->value);
}

int main() {
    fmd_table_t *table;
    fmd_table_init(&table, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    int no_items = 10000;
    for(int i = 0; i <= no_items; i++) {
        fmd_bs_t *bs;
        fmd_bs_init_reserve(&bs, 3);
        fmd_bs_init_write_word(bs,0,i*i,32);
        fmd_table_insert(table, i, bs);
    }
    fmd_table_rehash(&table, table->capacity*2);

    for(int i = no_items; i >= 0; i--) {
        fmd_bs_t *bs = 0;
        bool a = fmd_table_lookup(table, i, &bs);
        if(!a) exit(-1);
        word_t val;
        fmd_bs_init_read_word(bs,0,32, &val);
        printf("%d, %d\n", i, val);
    }
    return 0;
}