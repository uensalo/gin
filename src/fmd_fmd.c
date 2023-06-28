#include "fmd_fmd.h"
#include "fmd_tree.h"

fmd_fork_node_t *fmd_fork_node_init(fmd_fork_node_t *parent, int_t vlo, int_t vhi, int_t pos, bool is_leaf, bool is_dead, bool is_cadet) {
    fmd_fork_node_t *fn = calloc(1, sizeof(fmd_fork_node_t));
    fn->parent = (struct fmd_fork_node_t*)parent;
    fn->vertex_lo = vlo;
    fn->vertex_hi = vhi;
    fn->pos = pos;
    fn->is_leaf = is_leaf;
    fn->is_dead = is_dead;
    fn->is_cadet = is_cadet;
    return fn;
}

void fmd_fork_node_free(fmd_fork_node_t *node) {
    // frees only one node
    if(node) free(node);
}

fmd_fork_node_t *fmd_fork_node_copy(fmd_fork_node_t *node) {
    if(!node) return NULL;
    fmd_fork_node_t *copy = calloc(1, sizeof(fmd_fork_node_t));
    if(!copy) {
        free(copy);
        return NULL;
    }
    memcpy(copy, node, sizeof(fmd_fork_node_t));
    return copy;
}
uint_t fmd_fork_node_hash(fmd_fork_node_t *node) {
    const uint_t prime = 1099511628211LLU;
    uint_t hash = 14695981039346656037LLU;
    hash ^= prm_hash_f((void*)node->vertex_lo);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->vertex_hi);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->sa_lo);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->sa_hi);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->pos);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->is_leaf);
    hash *= prime;
    return hash;
}

int fmd_fork_node_comp(fmd_fork_node_t *n1, fmd_fork_node_t *n2) {
    return memcmp(n1, n2, sizeof(fmd_fork_node_t));
}

fmd_fmd_qr_t *fmd_fmd_qr_init(fmd_fork_node_t *cur_fork, int_t lo, int_t hi, int_t pos, fmd_string_t *pattern) {
    fmd_fmd_qr_t *qr = calloc(1, sizeof(fmd_fmd_qr_t));
    qr->cur_fork = cur_fork;
    qr->lo = lo;
    qr->hi = hi;
    qr->pos = pos;
    qr->pattern = pattern;
    return qr;
}

void fmd_fmd_qr_free(fmd_fmd_qr_t *q) {
    if(q) free(q);
}

void fmd_fmd_init(fmd_fmd_t** fmd, fmd_graph_t *graph, fmd_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate) {
    fmd_fmd_t *f = calloc(1, sizeof(fmd_fmd_t));
    if(!f || !graph) {
        *fmd = NULL;
        return;
    }
    f->c_0 = c_0;
    f->c_1 = c_1;
    /**************************************************************************
    * Step 0 - If no permutation is provided, generate an identity permutation
    **************************************************************************/
    if(!permutation) {
        fmd_vector_init(&permutation, graph->vertex_list->size, &prm_fstruct);
        for(int_t i = 0; i < graph->vertex_list->size; i++) {
            fmd_vector_append(permutation, (void*)i);
        }
    }
    assert(permutation->size == graph->vertex_list->size);
    assert(c_0 < c_1);
    f->permutation = permutation;
    /**************************************************************************
    * Step 1 - Compute the graph encoding and the permutation encodings
    **************************************************************************/
    /******************************************************
    * Step 1a - Generate fixed length binary encodings
    ******************************************************/
    fmd_vector_t *cwords = fmd_fmd_init_pcodes_fixed_binary_helper(FMD_FMD_DEFAULT_a_0, FMD_FMD_DEFAULT_a_1, graph->vertex_list->size);
    /******************************************************
    * Step 1b - Map the inverse of perm. to codewords
    ******************************************************/
    fmd_vector_t *inverse_permutation;
    fmd_vector_init(&inverse_permutation, permutation->size, &prm_fstruct);
    //fmd_table_t *cword_map;
    //fmd_table_init(&cword_map, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    // if(!cword_map...)
    for(int_t i = 0; i < permutation->size; i++) {
        inverse_permutation->data[(int_t)permutation->data[i]] = (void*)i;
        //fmd_table_insert(cword_map, permutation->data[i], (void*)i);
    }
    /******************************************************
    * Step 1c - Finally, create the graph encoding
    ******************************************************/
    // precompute size requirements
    int_t V = graph->vertex_list->size;
    int_t total_label_len = 0;
    for(int_t i = 0; i < V; i++) {
        total_label_len += ((fmd_vertex_t*)graph->vertex_list->data[i])->label->size;
    }
    // concatenate everything
    fmd_string_t *graph_encoding;
    fmd_string_init(&graph_encoding, total_label_len + (2 + fmd_ceil_log2(V)) * V);
    for(int_t i = 0; i < V; i++) {
        fmd_string_append(graph_encoding, c_0);
        fmd_string_concat_mut(graph_encoding, ((fmd_vertex_t*)graph->vertex_list->data[i])->label);
        fmd_string_append(graph_encoding, c_1);
        int_t pcode_idx = (int_t)inverse_permutation->data[i];
        //fmd_table_lookup(cword_map, (void*)i, (void*)&pcode_idx);
        fmd_string_concat_mut(graph_encoding, (fmd_string_t*)(cwords->data[pcode_idx]));
    }
    // free obsolete stuff
    fmd_vector_free(cwords);
    //fmd_table_free(cword_map);
    /**************************************************************************
    * Step 2 - Compute the suffix array of the encoding and build an FMI
    **************************************************************************/
    fmd_fmi_t *fmi;
    int64_t *sa = calloc(graph_encoding->size+1, sizeof(uint64_t));
    // if(!sa...)
    divsufsort64((sauchar_t*)graph_encoding->seq, (saidx64_t*)sa, graph_encoding->size+1);
    fmd_fmi_init_with_sa(&fmi,graph_encoding,sa,rank_sample_rate,isa_sample_rate);
    // if(!fmi...)
    f->graph_fmi = fmi;
    // DEBUG PURPOSES
    /*
    printf("S: %s%c\n", graph_encoding->seq,'\0');
    printf("F: ");
    for(int_t i = 0; i < graph_encoding->size+1; i++) {
        printf("%c", graph_encoding->seq[sa[i]]);
    }
    printf("\n");
    // DEBUG PURPOSES
    printf("L: ");
    for(int_t i = 0; i < graph_encoding->size+1; i++) {
        if (sa[i]) printf("%c", graph_encoding->seq[sa[i]-1]);
        else printf("%c", '\0');
    }
    printf("\n");
    printf("i:  ");
    for(int_t i = 0; i < graph_encoding->size+1; i++) {
        printf("%d ", i);
    }
    printf("\n");
    printf("SA: ");
    for(int_t i = 0; i < graph_encoding->size+1; i++) {
        printf("%d ", sa[i]);
    }
    printf("\n");
    */
    // END DEBUG PURPOSES

    /**************************************************************************
    * Step 3 - Compute the query range translation index
    **************************************************************************/
    /******************************************************
    * Step 3a - Get slices of char ranges of c_0, c_1
    ******************************************************/
    // sa[0] -> suffix starting with \0, essentially the terminator
    // sa[1]   ... sa[V]  : range of c_0
    // sa[V+1] ... sa[2V] : range of c_1
    fmd_vector_t *c_0_bucket, *c_1_bucket;
    fmd_vector_init(&c_0_bucket, V, &prm_fstruct);
    fmd_vector_init(&c_1_bucket, V, &prm_fstruct);
    // split in two loops for better cache performance
    for(int_t i = 1; i <= V; i++) {
        //c_0_bucket->data[i] = (void*)sa[i];
        fmd_vector_append(c_0_bucket, (void*)sa[i]);
        // DEBUG PURPOSES
        //assert(graph_encoding->seq[sa[i]] == c_0);
    }
    for(int_t i = V+1; i <= 2*V; i++) {
        fmd_vector_append(c_1_bucket, (void*)sa[i]);
        // DEBUG PURPOSES
        //assert(graph_encoding->seq[sa[i]] == c_1);
    }
    // at this point, sa is no longer needed as slices are extracted
    fmd_string_free(graph_encoding);
    free(sa); // frees gigs of memory :)
    /******************************************************
    * Step 3b - Compute rank translation tables
    ******************************************************/
    // get bwt ranks to text ranks mappings
    fmd_vector_t *c_0_bwt_to_text_tmp, *c_1_bwt_to_text, *c_1_text_to_bwt;
    fmd_vector_argsort(&c_0_bwt_to_text_tmp, c_0_bucket);
    fmd_vector_argsort(&c_1_bwt_to_text, c_1_bucket);


    // DEBUG PURPOSES
    fmd_vector_t *c_0_bwt_to_text;
    fmd_vector_argsort(&c_0_bwt_to_text, c_0_bwt_to_text_tmp);
    fmd_vector_free(c_0_bwt_to_text_tmp);
    // END OF DEBUG PURPOSES

    // now, invert the mapping of c_1 to get a mapping from text to bwt for c_1 ranks
    fmd_vector_init(&c_1_text_to_bwt, V, &prm_fstruct);
    fmd_table_t *c_1_inversion_table;
    fmd_table_init(&c_1_inversion_table, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    for(int_t i = 0; i < V; i++) {
        fmd_table_insert(c_1_inversion_table, (void*)c_1_bwt_to_text->data[i], (void*)i);
    }
    for(int_t i = 0; i < V; i++) {
        int_t val = -1;
        bool found = fmd_table_lookup(c_1_inversion_table, (void*)i, (void*)&val);
        assert(found);
        fmd_vector_append(c_1_text_to_bwt, (void*)val);
    }
    // DEBUG PURPOSES
    //for(int_t i = 0; i < c_1_text_to_bwt->size; i++) {
    //    assert(c_1_text_to_bwt->data[i] == permutation->data[i]);
    //}

    //printf("C_0 BWT to text:\n");
    //for(int_t i = 0; i < c_0_bwt_to_text->size; i++) {
    //    printf("%d:%d\n",i,c_0_bwt_to_text->data[i]);
    //}
    //printf("C_1 text to BWT:\n");
    //for(int_t i = 0; i < c_1_bwt_to_text->size; i++) {
    //    printf("%d:%d\n",i,c_1_text_to_bwt->data[i]);
    //}
    // END DEBUG PURPOSES

    // free intermediate structures
    fmd_vector_free(c_0_bucket);
    fmd_vector_free(c_1_bucket);
    fmd_vector_free(c_1_bwt_to_text);
    fmd_table_free(c_1_inversion_table);
    /******************************************************
    * Step 3c - (key,value)s for the interval merge tree
    ******************************************************/
    // keys : integers from 1 to V (vid+1)
    // values : list of intervals
    fmd_vector_t *kv_pairs;
    fmd_vector_init(&kv_pairs, V, &fmd_fstruct_vector); // vector holds interval lists
    for(int_t i = 0; i < V; i++) {
        int_t vid = (int_t)c_0_bwt_to_text->data[i];
        fmd_vector_t *neighbors;
        bool found = fmd_table_lookup(graph->incoming_neighbors, (void*)vid, (void*)&neighbors);
        assert(found);
        fmd_vector_t *neighbors_bwt_idx;
        fmd_vector_init(&neighbors_bwt_idx, neighbors->size, &prm_fstruct);
        for(int_t j = 0; j < neighbors->size; j++) {
            fmd_vector_append(neighbors_bwt_idx, (void*)inverse_permutation->data[(vid_t)neighbors->data[j]]);
        }
        fmd_vector_t *bwt_neighbor_intervals;
        fmd_vector_init(&bwt_neighbor_intervals, neighbors->size, &fmd_fstruct_imt_interval);
        if(neighbors_bwt_idx->size) {
            // DEBUG PURPOSES
            /*
            printf("bwt-range of %d: ", vid);
            for(int_t k = 0; k < neighbors_bwt_idx->size; k++) {
                printf("%d ", neighbors_bwt_idx->data[k]);
            }
            printf("\n");
             */
             // END OF DEBUG PURPOSES

            // sort and compact neighbors_bwt_idx into intervals
            fmd_vector_sort(neighbors_bwt_idx);
            // compaction start
            int_t lo = (int_t) neighbors_bwt_idx->data[0];
            int_t hi = lo;
            for (int_t j = 1; j < neighbors_bwt_idx->size; j++) {
                int_t current = (int_t) neighbors_bwt_idx->data[j];
                if (current == hi + 1) {
                    hi = current;
                } else {
                    fmd_imt_interval_t *interval;
                    fmd_imt_interval_init(&interval, lo, hi);
                    fmd_vector_append(bwt_neighbor_intervals, interval);
                    lo = hi = current;
                }
            }
            // don't forget to add the last interval
            fmd_imt_interval_t *interval;
            fmd_imt_interval_init(&interval, lo, hi);
            fmd_vector_append(bwt_neighbor_intervals, interval);
        }
        // free the singular list
        fmd_vector_free(neighbors_bwt_idx);
        fmd_vector_append(kv_pairs, bwt_neighbor_intervals);
    }
    f->bwt_to_vid = c_0_bwt_to_text;
    fmd_vector_free(c_1_text_to_bwt);
    /******************************************************
    * Step 3d - Now, actually construct the tree
    ******************************************************/
    // DEBUG PURPOSES
    /*
    for(int_t i = 0; i < kv_pairs->size; i++) {
        printf("%d : [ ", i);
        fmd_vector_t *intervals = kv_pairs->data[i];
        for(int_t j = 0; j < intervals->size; j++) {
            fmd_imt_interval_t *interval = intervals->data[j];
            printf("[%d,%d]  ",interval->lo, interval->hi);
        }
        printf("]\n");
    }
     */
    // END OF DEBUG PURPOSES
    fmd_imt_t *inverval_merge_tree;
    fmd_imt_init(&inverval_merge_tree, V, kv_pairs);
    f->r2r_tree = inverval_merge_tree;
    // if(!interval_merge_tree...
    // free the list structure storing interval lists, but not the lists themselves
    kv_pairs->f = &prm_fstruct;
    fmd_vector_free(kv_pairs);
    /******************************************************
    * Step 4 - Return
    ******************************************************/
    *fmd = f;
}

void fmd_fmd_free(fmd_fmd_t *fmd) {
    if(fmd) {
        fmd_vector_free(fmd->permutation);
        fmd_vector_free(fmd->bwt_to_vid);
        fmd_fmi_free(fmd->graph_fmi);
        fmd_imt_free(fmd->r2r_tree);
        free(fmd);
    }
}

int fmd_fmd_comp(fmd_fmd_t *f1, fmd_fmd_t *f2) {
    int c1 = fmd_fmi_comp(f1->graph_fmi, f2->graph_fmi);
    int c2 = fmd_imt_comp(f1->r2r_tree, f2->r2r_tree);
    int c3 = fmd_vector_comp(f1->permutation, f2->permutation);
    int c4 = fmd_vector_comp(f1->bwt_to_vid, f2->bwt_to_vid);
    int c5 = f1->c_0 == f2->c_0;
    int c6 = f1->c_1 == f2->c_1;
    return c1 == 0 && c2 == 0 && c3 == 0 && c4 == 0 && c5 && c6 ? 0 : -1;
}

count_t fmd_fmd_query_count(fmd_fmd_t *fmd, fmd_string_t *string) {
    // walk root count
    int_t count = 0;
    // push initial query
    fmd_vector_t *stack;
    fmd_vector_init(&stack, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    int_t V = fmd->permutation->size;
    int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fmi_qr_t *root_query = fmd_fmi_qr_init(init_lo, init_hi, string->size-1, string);
    fmd_vector_append(stack, root_query);
    while(stack->size > 0) {
        // pop one from the stack
        fmd_fmi_qr_t *query;
        fmd_vector_pop(stack, (void*)&query);
        while(query->pos > -1) {
            // now check: are there any exhausted vertices?
            int_t c_0_lo, c_0_hi;
            bool ok = fmd_fmi_query_precedence_range(fmd->graph_fmi, query, fmd->c_0, &c_0_lo, &c_0_hi);
            // if(!ok... // this check is not needed if everything goes ok in coding time :)
            // check if we need to fork into other vertices or not
            if(c_0_hi > c_0_lo) {
                // we have a walk having the current suffix of the query as a prefix
                fmd_vector_t *incoming_sa_intervals;
                fmd_imt_query(fmd->r2r_tree, c_0_lo-1, c_0_hi-2, &incoming_sa_intervals);
                for(int_t i = 0; i < incoming_sa_intervals->size; i++) {
                    fmd_imt_interval_t *interval = incoming_sa_intervals->data[i];
                    fmd_fmi_qr_t *fork = fmd_fmi_qr_init(V+1+interval->lo, V+interval->hi+2, query->pos, string);
                    fmd_vector_append(stack, fork);
                }
                fmd_vector_free(incoming_sa_intervals);
            }
            bool okc = fmd_fmi_advance_query(fmd->graph_fmi, query);
            if(!okc || query->lo == query->hi) {
                break;
            }
        }
        count += query->hi - query->lo;
        fmd_fmi_qr_free(query);
    }
    fmd_vector_free(stack);
    return count;
}

fmd_vector_t *fmd_fmd_query_locate_basic(fmd_fmd_t *fmd, fmd_string_t *string) {
    // walk root count
    fmd_vector_t *locs;
    fmd_vector_init(&locs, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    // push initial query
    fmd_vector_t *stack;
    fmd_vector_init(&stack, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    int_t V = fmd->permutation->size;
    int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fmi_qr_t *root_query = fmd_fmi_qr_init(init_lo, init_hi, string->size-1, string);
    fmd_vector_append(stack, root_query);
    while(stack->size > 0) {
        // pop one from the stack
        fmd_fmi_qr_t *query;
        fmd_vector_pop(stack, (void*)&query);
        while(query->pos > -1) {
            // now check: are there any exhausted vertices?
            int_t c_0_lo, c_0_hi;
            bool okc = fmd_fmi_advance_query(fmd->graph_fmi, query);
            bool ok = fmd_fmi_query_precedence_range(fmd->graph_fmi, query, fmd->c_0, &c_0_lo, &c_0_hi);
            // if(!ok... // this check is not needed if everything goes ok in coding time :)
            if (query->pos == -1) break;
            // check if we need to fork into other vertices or not
            if(c_0_hi > c_0_lo) {
                // we have a walk having the current suffix of the query as a prefix
                fmd_vector_t *incoming_sa_intervals;
                fmd_imt_query(fmd->r2r_tree, c_0_lo-1, c_0_hi-2, &incoming_sa_intervals);
                for(int_t i = 0; i < incoming_sa_intervals->size; i++) {
                    fmd_imt_interval_t *interval = incoming_sa_intervals->data[i];
                    fmd_fmi_qr_t *fork = fmd_fmi_qr_init(V+1+interval->lo, V+interval->hi+2, query->pos, string);
                    fmd_vector_append(stack, fork);
                }
                fmd_vector_free(incoming_sa_intervals);
            }
            if(!okc || query->lo == query->hi) {
                break;
            }
        }
        fmd_vector_t *qlocs = fmd_fmi_sa(fmd->graph_fmi, query);
        for(int_t i = 0; i < qlocs->size; i++) {
            fmd_vector_append(locs, qlocs->data[i]);
        }
        fmd_vector_free(qlocs);
        fmd_fmi_qr_free(query);
    }
    fmd_vector_free(stack);
    return locs;
}

void fmd_fmd_query_locate_paths(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **paths, fmd_vector_t **dead_ends) {
    // walk root count
    fmd_vector_t *leaves;
    fmd_vector_init(&leaves, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    fmd_vector_t *graveyard;
    fmd_vector_init(&graveyard, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    // push initial query
    fmd_vector_t *stack;
    fmd_vector_init(&stack, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    int_t V = fmd->permutation->size;
    int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fork_node_t *root_fork = fmd_fork_node_init(NULL,
                                                    -1, -1,
                                                    string->size-1,
                                                    false, false, false);
    fmd_fmd_qr_t *root_query = fmd_fmd_qr_init(root_fork, init_lo, init_hi, string->size-1, string);
    fmd_vector_append(stack, root_query);
    while(stack->size > 0) {
        // pop one from the stack
        fmd_fmd_qr_t *query;
        fmd_vector_pop(stack, (void*)&query);
        while(query->pos > -1) {
            // now check: are there any exhausted vertices?
            int_t c_0_lo, c_0_hi;
            bool okc = fmd_fmd_advance_query(fmd->graph_fmi, query);
            bool ok = fmd_fmd_query_precedence_range(fmd->graph_fmi, query, fmd->c_0, &c_0_lo, &c_0_hi);
            // if(!ok... // this check is not needed if everything goes ok in coding time :)
            if (query->pos == -1) break;
            // check if we need to fork into other vertices or not
            if(c_0_hi > c_0_lo) {
                // we have a walk having the current suffix of the query as a prefix
                fmd_fork_node_t *royal_node = fmd_fork_node_init(query->cur_fork, query->cur_fork->vertex_lo, query->cur_fork->vertex_hi, query->pos, false, false, false);
                fmd_vector_t *incoming_sa_intervals;
                fmd_imt_query(fmd->r2r_tree, c_0_lo-1, c_0_hi-2, &incoming_sa_intervals);
                for(int_t i = 0; i < incoming_sa_intervals->size; i++) {
                    fmd_imt_interval_t *interval = incoming_sa_intervals->data[i];
                    fmd_fork_node_t *cadet_node = fmd_fork_node_init(query->cur_fork, interval->lo, interval->hi+1, query->pos, false, false, true);
                    fmd_fmd_qr_t *fork = fmd_fmd_qr_init(cadet_node, V+1+interval->lo, V+2+interval->hi, query->pos, string);
                    fmd_vector_append(stack, fork);
                }
                query->cur_fork = royal_node;
                fmd_vector_free(incoming_sa_intervals);
            }
            if(!okc || query->lo == query->hi) {
                fmd_fork_node_t *dead_node = fmd_fork_node_init(query->cur_fork, c_0_lo, c_0_hi, query->pos, true, true, false);
                dead_node->sa_lo = query->lo;
                dead_node->sa_hi = query->hi;
                query->cur_fork = dead_node;
                break;
            }
        }
        if(!query->cur_fork->is_dead && !query->cur_fork->is_leaf) {
            fmd_fork_node_t *leaf_node = fmd_fork_node_init(query->cur_fork, query->cur_fork->vertex_lo,
                                                            query->cur_fork->vertex_hi, query->pos, true, false, false);
            leaf_node->sa_lo = query->lo;
            leaf_node->sa_hi = query->hi;
            if(query->lo >= query->hi) {
                leaf_node->is_dead = true;
                fmd_vector_append(graveyard, leaf_node);
            } else {
                fmd_vector_append(leaves, leaf_node);
            }
        } else {
            fmd_vector_append(graveyard, query->cur_fork);
        }
        fmd_fmd_qr_free(query);
    }
    fmd_vector_free(stack);
    *paths = leaves;
    *dead_ends = graveyard;
}

void fmd_fmd_locate_paths_result_free(fmd_vector_t *paths, fmd_vector_t *dead_ends) {
    // keys are pointers, values are actual pointers
    fmd_tree_t *visited;
    fmd_tree_init(&visited, &prm_fstruct, &fmd_fstruct_fork_node);

    for(int_t i = 0; i < paths->size; i++) {
        fmd_fork_node_t *cur = paths->data[i];
        while(cur) {
            bool inserted = fmd_tree_insert(visited, cur, cur);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    for(int_t i = 0; i < dead_ends->size; i++) {
        fmd_fork_node_t *cur = dead_ends->data[i];
        while(cur) {
            bool inserted = fmd_tree_insert(visited, cur, cur);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    fmd_tree_free(visited);
    paths->f = &prm_fstruct;
    fmd_vector_free(paths);
    dead_ends->f = &prm_fstruct;
    fmd_vector_free(dead_ends);
}

void fmd_fmd_query_locate_paths_stats(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **paths, fmd_vector_t **dead_ends, int_t *no_forks) {
    int_t forks = 0;
    // walk root count
    fmd_vector_t *leaves;
    fmd_vector_init(&leaves, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    fmd_vector_t *graveyard;
    fmd_vector_init(&graveyard, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    // push initial query
    fmd_vector_t *stack;
    fmd_vector_init(&stack, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    int_t V = fmd->permutation->size;
    int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fork_node_t *root_fork = fmd_fork_node_init(NULL,
                                                    -1, -1,
                                                    string->size-1,
                                                    false, false, false);
    fmd_fmd_qr_t *root_query = fmd_fmd_qr_init(root_fork, init_lo, init_hi, string->size-1, string);
    fmd_vector_append(stack, root_query);
    while(stack->size > 0) {
        // pop one from the stack
        fmd_fmd_qr_t *query;
        fmd_vector_pop(stack, (void*)&query);
        while(query->pos > -1) {
            // now check: are there any exhausted vertices?
            int_t c_0_lo, c_0_hi;
            bool okc = fmd_fmd_advance_query(fmd->graph_fmi, query);
            bool ok = fmd_fmd_query_precedence_range(fmd->graph_fmi, query, fmd->c_0, &c_0_lo, &c_0_hi);
            // if(!ok... // this check is not needed if everything goes ok in coding time :)
            if (query->pos == -1) break;
            // check if we need to fork into other vertices or not
            if(c_0_hi > c_0_lo) {
                // we have a walk having the current suffix of the query as a prefix
                fmd_fork_node_t *royal_node = fmd_fork_node_init(query->cur_fork, query->cur_fork->vertex_lo, query->cur_fork->vertex_hi, query->pos, false, false, false);
                fmd_vector_t *incoming_sa_intervals;
                fmd_imt_query(fmd->r2r_tree, c_0_lo-1, c_0_hi-2, &incoming_sa_intervals);
                for(int_t i = 0; i < incoming_sa_intervals->size; i++) {
                    fmd_imt_interval_t *interval = incoming_sa_intervals->data[i];
                    fmd_fork_node_t *cadet_node = fmd_fork_node_init(query->cur_fork, interval->lo, interval->hi+1, query->pos, false, false, true);
                    fmd_fmd_qr_t *fork = fmd_fmd_qr_init(cadet_node, V+1+interval->lo, V+2+interval->hi, query->pos, string);
                    fmd_vector_append(stack, fork);
                }
                forks += incoming_sa_intervals->size;
                query->cur_fork = royal_node;
                fmd_vector_free(incoming_sa_intervals);
            }
            if(!okc || query->lo == query->hi) {
                fmd_fork_node_t *dead_node = fmd_fork_node_init(query->cur_fork, c_0_lo, c_0_hi, query->pos, true, true, false);
                dead_node->sa_lo = query->lo;
                dead_node->sa_hi = query->hi;
                query->cur_fork = dead_node;
                break;
            }
        }
        if(!query->cur_fork->is_dead && !query->cur_fork->is_leaf) {
            fmd_fork_node_t *leaf_node = fmd_fork_node_init(query->cur_fork, query->cur_fork->vertex_lo,
                                                            query->cur_fork->vertex_hi, query->pos, true, false, false);
            leaf_node->sa_lo = query->lo;
            leaf_node->sa_hi = query->hi;
            if(query->lo >= query->hi) {
                leaf_node->is_dead = true;
                fmd_vector_append(graveyard, leaf_node);
            } else {
                fmd_vector_append(leaves, leaf_node);
            }
        } else {
            fmd_vector_append(graveyard, query->cur_fork);
        }
        fmd_fmd_qr_free(query);
    }
    fmd_vector_free(stack);
    *paths = leaves;
    *dead_ends = graveyard;
    *no_forks = forks;
}

bool fmd_fmd_advance_query(fmd_fmi_t *fmi, fmd_fmd_qr_t *qr) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    bool found = fmd_table_lookup(fmi->c2e, (void*)qr->pattern->seq[qr->pos], &encoding);
    if(!found) {
        qr->lo = 0;
        qr->hi = 0;
        //fprintf(stderr,"[fmd_fmi_advance_query]: encoding not found in dictionary, query is NIL\n");
        return false;
    }
    count_t rank_lo_m_1 = qr->lo ? fmd_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? fmd_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    qr->lo = (int_t)(base + rank_lo_m_1);
    qr->hi = (int_t)(base + rank_hi_m_1);
    --qr->pos;
    return true;
}

bool fmd_fmd_query_precedence_range(fmd_fmi_t *fmi, fmd_fmd_qr_t *qr, char_t c, int_t *lo, int_t *hi) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    bool found = fmd_table_lookup(fmi->c2e, (void*)c, &encoding);
    if(!found) {
        qr->lo = 0;
        qr->hi = 0;
        //fprintf(stderr,"[fmd_fmi_advance_query]: encoding not found in dictionary, query is NIL\n");
        return false;
    }
    count_t rank_lo_m_1 = qr->lo ? fmd_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? fmd_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    *lo = (int_t)(base + rank_lo_m_1);
    *hi = (int_t)(base + rank_hi_m_1);
    return true;
}

fmd_vector_t *fmd_fmd_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords) {
    fmd_vector_t *rval;
    fmd_vector_init(&rval,no_codewords,&fmd_fstruct_string);
    int_t len = (int_t)fmd_ceil_log2(no_codewords); // length of each codeword
    char codeword[len+1]; // buffer for each generated codeword (plus a null terminator)
    for (int_t i = 0; i < no_codewords; i++) {
        int_t num = i; // binary representation of the current number
        int_t pos = len-1; // current position in the codeword buffer
        while (pos >= 0) {
            if (num % 2 == 0) {
                codeword[pos] = a_0;
            } else {
                codeword[pos] = a_1;
            }
            num /= 2;
            pos--;
        }
        codeword[len] = '\0'; // null terminate the codeword string
        fmd_string_t *copy;
        fmd_string_init_cstr(&copy, codeword);
        fmd_vector_append(rval, copy);
    }
    return rval;
}

void fmd_fmd_serialize_from_buffer(fmd_fmd_t **fmd_ret, unsigned char *buf, uint64_t buf_size) {
    fmd_fmd_t *fmd = calloc(1, sizeof(fmd_fmd_t));
    if(!fmd) {
        *fmd_ret = NULL;
        return;
    }
    /**************************************************************************
    * Step 0 - Wrap the buffer into a bitstream - this consumes the buffer
    **************************************************************************/
    fmd_bs_t *bs;
    fmd_bs_init_from_buffer(buf, buf_size, &bs);
    uint_t ridx = 0;
    /**************************************************************************
    * Step 1 - Read total size and special characters
    **************************************************************************/
    word_t fmd_size_in_bits, c_0, c_1;
    fmd_bs_read_word(bs, ridx, FMD_FMD_NO_BITS_BIT_LENGTH, &fmd_size_in_bits);
    ridx += FMD_FMD_NO_BITS_BIT_LENGTH;
    fmd_bs_read_word(bs, ridx, FMD_FMD_SPECIAL_CHAR_BIT_LENGTH, &c_0);
    ridx += FMD_FMD_SPECIAL_CHAR_BIT_LENGTH;
    fmd_bs_read_word(bs, ridx, FMD_FMD_SPECIAL_CHAR_BIT_LENGTH, &c_1);
    ridx += FMD_FMD_SPECIAL_CHAR_BIT_LENGTH;
    fmd->c_0 = (char_t)c_0;
    fmd->c_1 = (char_t)c_1;
    /**************************************************************************
    * Step 2 - Read the permutation and BWT to vertex ID map
    **************************************************************************/
    word_t no_vertices;
    fmd_bs_read_word(bs, ridx, FMD_FMD_NO_VERTICES_BIT_LENGTH, &no_vertices);
    ridx += FMD_FMD_NO_VERTICES_BIT_LENGTH;

    fmd_vector_t *permutation, *bwt_to_vid;
    fmd_vector_init(&permutation, (int_t)no_vertices, &prm_fstruct);
    fmd_vector_init(&bwt_to_vid, (int_t)no_vertices, &prm_fstruct);
    word_t p;
    for(int_t i = 0; i < no_vertices; i++) {
        fmd_bs_read_word(bs, ridx, FMD_FMD_PERMUTATION_BIT_LENGTH, &p);
        ridx += FMD_FMD_PERMUTATION_BIT_LENGTH;
        fmd_vector_append(permutation, (void*)p);
    }
    for(int_t i = 0; i < no_vertices; i++) {
        fmd_bs_read_word(bs, ridx, FMD_FMD_BWT_TO_VID_BIT_LENGTH, &p);
        ridx += FMD_FMD_BWT_TO_VID_BIT_LENGTH;
        fmd_vector_append(bwt_to_vid, (void*)p);
    }
    fmd_vector_fit(permutation);
    fmd_vector_fit(bwt_to_vid);
    fmd->permutation = permutation;
    fmd->bwt_to_vid = bwt_to_vid;
    /**************************************************************************
    * Step 3 - Reconstruct the FMI from the bitstream
    **************************************************************************/
    word_t fmi_no_bits;
    fmd_bs_read_word(bs, ridx, FMD_FMD_FMI_NO_BITS_BIT_LENGTH, &fmi_no_bits);
    ridx += FMD_FMD_FMI_NO_BITS_BIT_LENGTH;
    fmd_fmi_t *graph_fmi = calloc(1, sizeof(fmd_fmi_t));
    uint_t fmi_ridx_start = ridx;
    /******************************************************
    * Step 3a - Read the header
    ******************************************************/
    word_t no_chars, no_chars_per_block, isa_sample_rate, alphabet_size;
    fmd_bs_read_word(bs, ridx, FMD_FMI_CHAR_COUNT_BIT_LENGTH, &no_chars);
    ridx += FMD_FMI_CHAR_COUNT_BIT_LENGTH;

    fmd_bs_read_word(bs, ridx, FMD_FMI_RNK_SAMPLE_RATE_BIT_LENGTH, &no_chars_per_block);
    ridx += FMD_FMI_RNK_SAMPLE_RATE_BIT_LENGTH;

    fmd_bs_read_word(bs, ridx, FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH, &isa_sample_rate);
    ridx += FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;

    fmd_bs_read_word(bs, ridx, FMD_FMI_ALPHABET_SIZE_BIT_LENGTH, &alphabet_size);
    ridx += FMD_FMI_ALPHABET_SIZE_BIT_LENGTH;

    graph_fmi->no_chars = (int_t)no_chars;
    graph_fmi->no_chars_per_block = (int_t)no_chars_per_block;
    graph_fmi->isa_sample_rate = (int_t)isa_sample_rate;
    graph_fmi->alphabet_size = (int_t)alphabet_size;
    graph_fmi->no_bits_per_char = fmd_ceil_log2((int_t)alphabet_size);
    /******************************************************
    * Step 3b - Read the alphabet
    ******************************************************/
    fmd_vector_init(&graph_fmi->alphabet, graph_fmi->alphabet_size, &prm_fstruct);
    fmd_table_init(&graph_fmi->e2c, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    fmd_table_init(&graph_fmi->c2e, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    for(int_t i = 0; i < graph_fmi->alphabet_size; i++) {
        word_t alphabet_char, encoding;
        fmd_bs_read_word(bs, ridx, FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH, &alphabet_char);
        ridx+=FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH;
        fmd_vector_append(graph_fmi->alphabet, (void*)alphabet_char);

        fmd_bs_read_word(bs, ridx, FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH, &encoding);
        ridx+=FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH;

        fmd_table_insert(graph_fmi->c2e, (void*)alphabet_char, (void*)encoding);
        fmd_table_insert(graph_fmi->e2c, (void*)encoding, (void*)alphabet_char);
    }
    /******************************************************
    * Step 3c - Read the suffix array samples
    ******************************************************/
    fmd_table_init(&graph_fmi->isa, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    word_t sa_idx;
    for(int_t i = 0; i < graph_fmi->no_chars / graph_fmi->isa_sample_rate; i++) {
        int_t sa_val = i * graph_fmi->isa_sample_rate;
        fmd_bs_read_word(bs, ridx, FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH, &sa_idx);
        ridx+=FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;
        fmd_table_insert(graph_fmi->isa, (void*)sa_idx, (void*)sa_val);
    }
    /******************************************************
    * Step 3d - Copy the bits field and the rest
    ******************************************************/
    graph_fmi->bv_start_offset = (int_t)ridx - (int_t)fmi_ridx_start;
    fmd_bs_t *fmi_bits;
    fmd_bs_init(&fmi_bits);
    fmd_bs_fit(fmi_bits, fmi_no_bits);
    for(int_t i = 0; i < fmi_bits->cap_in_words; i++) {
        word_t read_bits;
        uint_t offset = i * WORD_NUM_BITS;
        fmd_bs_read_word(bs, fmi_ridx_start + offset, WORD_NUM_BITS, &read_bits);
        fmd_bs_write_word(fmi_bits, offset, read_bits, WORD_NUM_BITS);
    }
    fmd_bs_fit(fmi_bits, fmi_no_bits);
    graph_fmi->bits = fmi_bits;
    graph_fmi->no_bits = (int_t)fmi_no_bits;
    /******************************************************
    * Step 3e - Compute the cumulative sum from last block
    ******************************************************/
    graph_fmi->char_counts = calloc(graph_fmi->alphabet_size, sizeof(count_t));
    count_t cum = 0;
    //uint_t block_idx = graph_fmi->no_chars / graph_fmi->no_chars_per_block;
    //uint_t block_size = (graph_fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + graph_fmi->no_bits_per_char * graph_fmi->no_chars_per_block);
    uint_t rank_cache_idx = graph_fmi->no_bits - graph_fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    for(int_t i = 0; i < graph_fmi->alphabet_size; i++) {
        word_t count;
        graph_fmi->char_counts[i] = cum;
        fmd_bs_read_word(graph_fmi->bits, rank_cache_idx + i * FMD_FMI_CHAR_COUNT_BIT_LENGTH, FMD_FMI_CHAR_COUNT_BIT_LENGTH, &count);
        cum += (count_t)count;
    }
    /******************************************************
    * Step 3f - Finalize reading the FMI
    ******************************************************/
    /* Don't forget to set ridx accordingly */
    fmd->graph_fmi = graph_fmi;
    ridx = fmi_ridx_start + fmi_bits->cap_in_words * WORD_NUM_BITS;

    /**************************************************************************
    * Step 4 - Construct the interval-merge tree from the leaves
    **************************************************************************/
    word_t fmd_imt_bit_length;
    fmd_bs_read_word(bs, ridx, FMD_FMD_IMT_NO_BITS_BIT_LENGTH, &fmd_imt_bit_length);
    ridx += FMD_FMD_IMT_NO_BITS_BIT_LENGTH;
    // no keys is equal to the number of vertices
    fmd_vector_t *kv_pairs;
    fmd_vector_init(&kv_pairs, (int_t)no_vertices, &fmd_fstruct_vector);
    word_t no_intervals, lo, hi;
    for(int_t i = 0; i < no_vertices; i++) {
        fmd_bs_read_word(bs, ridx, FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH, &no_intervals);
        ridx += FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH;
        fmd_vector_t *bwt_neighbor_intervals;
        fmd_vector_init(&bwt_neighbor_intervals, (int_t)no_intervals, &fmd_fstruct_imt_interval);
        for(int_t j = 0; j < (int_t)no_intervals; j++) {
            fmd_bs_read_word(bs, ridx, FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH, &lo);
            ridx += FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            fmd_bs_read_word(bs, ridx, FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH, &hi);
            ridx += FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            fmd_imt_interval_t *interval;
            fmd_imt_interval_init(&interval, (int_t)lo, (int_t)hi);
            fmd_vector_append(bwt_neighbor_intervals, interval);
        }
        fmd_vector_append(kv_pairs, bwt_neighbor_intervals);
    }
    fmd_imt_t *inverval_merge_tree;
    fmd_imt_init(&inverval_merge_tree, (int_t)no_vertices, kv_pairs);
    fmd->r2r_tree = inverval_merge_tree;
    kv_pairs->f = &prm_fstruct;
    fmd_vector_free(kv_pairs);
    /**************************************************************************
    * Step 5 - Cleanup and return the reconstructed index
    **************************************************************************/
    fmd_bs_free(bs);
    *fmd_ret = fmd;
}

void fmd_fmd_serialize_to_buffer(fmd_fmd_t *fmd, unsigned char **buf_ret, uint64_t *buf_size_ret) {
    if(!fmd) {
        *buf_ret = NULL;
        *buf_size_ret = 0;
        return;
    }
    // init bitstream to write everything
    fmd_bs_t *bs;
    fmd_bs_init(&bs);
    uint_t widx = 0;
    /**************************************************************************
    * Step 0 - Write the special characters to the buffer
    **************************************************************************/
    /* Reserve the beginning of the buffer for the total bit length */
    widx += FMD_FMD_NO_BITS_BIT_LENGTH;
    fmd_bs_write_word(bs, widx, fmd->c_0, FMD_FMD_SPECIAL_CHAR_BIT_LENGTH);
    widx += FMD_FMD_SPECIAL_CHAR_BIT_LENGTH;
    fmd_bs_write_word(bs, widx, fmd->c_1, FMD_FMD_SPECIAL_CHAR_BIT_LENGTH);
    widx += FMD_FMD_SPECIAL_CHAR_BIT_LENGTH;
    /**************************************************************************
    * Step 1 - Write the permutation to the buffer
    **************************************************************************/
    fmd_bs_write_word(bs, widx, fmd->permutation->size, FMD_FMD_NO_VERTICES_BIT_LENGTH);
    widx += FMD_FMD_NO_VERTICES_BIT_LENGTH;
    for(int_t i = 0; i < fmd->permutation->size; i++) {
        fmd_bs_write_word(bs, widx, (word_t)fmd->permutation->data[i], FMD_FMD_PERMUTATION_BIT_LENGTH);
        widx += FMD_FMD_PERMUTATION_BIT_LENGTH;
    }
    /**************************************************************************
    * Step 2 - Write the BWT c_0 to VID mapping
    **************************************************************************/
    for(int_t i = 0; i < fmd->bwt_to_vid->size; i++) {
        fmd_bs_write_word(bs, widx, (word_t)fmd->bwt_to_vid->data[i], FMD_FMD_BWT_TO_VID_BIT_LENGTH);
        widx += FMD_FMD_BWT_TO_VID_BIT_LENGTH;
    }
    /**************************************************************************
    * Step 3 - Write the FMI of the graph encoding
    **************************************************************************/
    /**************************************************************************
    * Step 4 - Write the FMI bitstream
    **************************************************************************/
    fmd_bs_write_word(bs, widx, (word_t)fmd->graph_fmi->no_bits, FMD_FMD_FMI_NO_BITS_BIT_LENGTH);
    widx += FMD_FMD_FMI_NO_BITS_BIT_LENGTH;
    fmd_bs_fit(fmd->graph_fmi->bits, fmd->graph_fmi->no_bits);
    for(int_t i = 0; i < fmd->graph_fmi->bits->cap_in_words; i++) {
        fmd_bs_write_word(bs, widx, (word_t)fmd->graph_fmi->bits->words[i], WORD_NUM_BITS);
        widx += WORD_NUM_BITS;
    }
    /**************************************************************************
    * Step 5 - Write the IMT bitstream
    **************************************************************************/
    /* Reserve space for the IMT encoding bit length */
    widx += FMD_FMD_IMT_NO_BITS_BIT_LENGTH;
    uint_t prev_widx = widx;
    fmd_fmd_serialize_to_buffer_imt_helper(fmd->r2r_tree->root, bs, &widx);
    fmd_bs_write_word(bs, prev_widx - FMD_FMD_IMT_NO_BITS_BIT_LENGTH, widx - prev_widx, FMD_FMD_IMT_NO_BITS_BIT_LENGTH);
    /**************************************************************************
    * Step 6 - Write total length in the beginning
    **************************************************************************/
    fmd_bs_write_word(bs, 0, widx, FMD_FMD_NO_BITS_BIT_LENGTH);
    /**************************************************************************
    * Step 7 - Return
    **************************************************************************/
    fmd_bs_fit(bs, widx);
    word_t *buf;
    uint_t no_words;
    fmd_bs_detach(bs, &buf, &no_words);
    *buf_ret = (unsigned char*)buf;
    *buf_size_ret = no_words * sizeof(word_t);
}

void fmd_fmd_serialize_to_buffer_imt_helper(fmd_imt_node_t *node, fmd_bs_t *bs, uint_t *widx) {
    if(node) {
        fmd_fmd_serialize_to_buffer_imt_helper((fmd_imt_node_t*)node->left, bs, widx);
        if (node->hi == node->lo) {
            fmd_bs_write_word(bs, *widx, (word_t)node->intervals->size, FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH);
            *widx += FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH;
            for(int_t i = 0; i < node->intervals->size; i++) {
                int_t lo = ((fmd_imt_interval_t*)(node->intervals->data[i]))->lo;
                int_t hi = ((fmd_imt_interval_t*)(node->intervals->data[i]))->hi;
                fmd_bs_write_word(bs, *widx, lo, FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH);
                *widx += FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
                fmd_bs_write_word(bs, *widx, hi, FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH);
                *widx += FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            }
        }
        fmd_fmd_serialize_to_buffer_imt_helper((fmd_imt_node_t*)node->right, bs, widx);
    }
}