#include "fmd_fmd.h"


fmd_fork_node_t *fmd_fork_node_init(fmd_fork_node_t *parent,
                                    int_t sa_lo, int_t sa_hi,
                                    int_t pos,
                                    fmd_fork_node_type_t type) {
    fmd_fork_node_t *fn = calloc(1, sizeof(fmd_fork_node_t));
    fn->parent = (struct fmd_fork_node_t*)parent;
    fn->sa_lo = sa_lo;
    fn->sa_hi = sa_hi;
    fn->pos = pos;
    fn->type = type;
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
    hash ^= prm_hash_f((void*)node->sa_lo);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->sa_hi);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->pos);
    hash *= prime;
    hash ^= prm_hash_f((void*)node->type);
    hash *= prime;
    return hash;
}

int fmd_fork_node_comp_exact(fmd_fork_node_t *n1, fmd_fork_node_t *n2) {
    return memcmp(n1, n2, sizeof(fmd_fork_node_t));
}
int fmd_fork_node_comp(fmd_fork_node_t *n1, fmd_fork_node_t *n2) {
    return n1->sa_lo < n2->sa_lo ? -1 : (n1->sa_lo > n2->sa_lo);
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
    f->permutation = fmd_vector_copy(permutation);
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
    fmd_vector_free(inverse_permutation);
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

#ifdef FMD_OMP
#include <omp.h>
#endif
void fmd_fmd_query_find_dfs(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends, int_t num_threads) {
    fmd_vector_t *leaves;
    fmd_vector_init(&leaves, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    fmd_vector_t *graveyard = NULL;
    fmd_vector_init(&graveyard, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);

    // create initial task to fire
    int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fork_node_t *root_fork = fmd_fork_node_init(NULL,
                                                    init_lo, init_hi,
                                                    string->size-1,
                                                    ROOT);
#ifdef FMD_OMP
    omp_set_num_threads((int)num_threads);
#endif
    #pragma omp parallel default(none) shared(fmd,root_fork,leaves,graveyard,max_forks,string)
    {
        #pragma omp single nowait
        {
            fmd_fmd_query_find_dfs_process_fork(fmd,root_fork,max_forks,string,leaves,graveyard);
        }
    }
    #pragma omp taskwait
    *paths = leaves;
    *dead_ends = graveyard;
}

void fmd_fmd_query_find_dfs_process_fork(fmd_fmd_t *fmd, fmd_fork_node_t *fork, int_t max_forks, fmd_string_t *pattern, fmd_vector_t *exact_matches, fmd_vector_t *partial_matches) {
    bool continue_flag = true;
    #pragma omp critical(exact_match)
    {
        if(max_forks != -1 && exact_matches->size >= max_forks) {
            continue_flag = false;
        }
    }
    int_t V = fmd->permutation->size;
    if (!continue_flag) {
        fmd_fork_node_free(fork);
        return;
    }
    while (fork->pos > -1) {
        int_t c_0_lo, c_0_hi;
        bool okc = fmd_fmd_advance_fork(fmd, fork, pattern);
        bool ok = fmd_fmd_fork_precedence_range(fmd, fork, fmd->c_0, &c_0_lo, &c_0_hi);
        if (fork->pos == -1) break;
        if (c_0_hi > c_0_lo) {
            // we have a walk having the current suffix of the query as a prefix
            fmd_fork_node_t *royal_node = fmd_fork_node_init(fork,
                                                             fork->sa_lo, fork->sa_hi,
                                                             //fork->v_lo, fork->v_hi,
                                                             fork->pos,
                                                             MAIN);
            fmd_vector_t *incoming_sa_intervals;
            fmd_imt_query(fmd->r2r_tree, c_0_lo - 1, c_0_hi - 2, &incoming_sa_intervals);
            for (int_t i = 0; i < incoming_sa_intervals->size; i++) {
                fmd_imt_interval_t *interval = incoming_sa_intervals->data[i];
                fmd_fork_node_t *new_fork = fmd_fork_node_init(royal_node, V+1+interval->lo, V+2+interval->hi,
                                                                 //interval->lo, interval->hi + 1,
                                                                 fork->pos, FALT);
                // fire subtask
                #pragma omp task default(none) shared(fmd, new_fork, exact_matches, partial_matches, max_forks, pattern)
                fmd_fmd_query_find_dfs_process_fork(fmd, new_fork, max_forks, pattern, exact_matches, partial_matches);
            }
            fork = royal_node;
            fmd_vector_free(incoming_sa_intervals);
        }
        if (!okc || fork->sa_lo == fork->sa_hi) {
            fmd_fork_node_t *dead_node = fmd_fork_node_init(fork, fork->sa_lo, fork->sa_hi,
                                                            //c_0_lo, c_0_hi,
                                                            fork->pos,
                                                            DEAD);
            fork = dead_node;
            break;
        }
    }
    if (fork->type != DEAD && fork->type != LEAF) {
        fmd_fork_node_t *leaf_node = fmd_fork_node_init(fork, fork->sa_lo, fork->sa_hi,
                                                        //fork->v_lo, fork->v_hi,
                                                        fork->pos,
                                                        LEAF);
        if (fork->sa_lo >= fork->sa_hi) {
            leaf_node->type = DEAD;
            #pragma omp critical(partial_match)
            {
                fmd_vector_append(partial_matches, leaf_node);
            }
        } else {
            #pragma omp critical(exact_match)
            {
                if(max_forks == -1 || exact_matches->size < max_forks) {
                    fmd_vector_append(exact_matches, leaf_node);
                } else {
                    fmd_vector_append(partial_matches, leaf_node);
                }
            }
        }
    } else {
        #pragma omp critical(partial_match)
        {
            fmd_vector_append(partial_matches, fork);
        }
    }
    //fmd_fmd_qr_free(query);
}

void fmd_fmd_query_find_step(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_forks, int_t *t, fmd_vector_t **cur_forks, fmd_vector_t **partial_matches) {
    fmd_vector_t *forks= *cur_forks;
    int_t V = fmd->permutation->size;
    /**********************************************************************
    * Step 1 - Forking phase: Fork each query at its current position
    **********************************************************************/
    fmd_vector_t *new_forks;
    fmd_vector_init(&new_forks, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
#pragma omp parallel for default(none) shared(forks, fmd, max_forks, new_forks, V)
    for (int_t i = 0; i < forks->size; i++) {
        fmd_fork_node_t *fork = forks->data[i];
        int_t c_0_lo, c_0_hi;
        fmd_fmd_fork_precedence_range(fmd, fork, fmd->c_0, &c_0_lo, &c_0_hi);
        bool more_to_track = max_forks == -1 || max_forks > forks->size;
        if(more_to_track && c_0_lo < c_0_hi) {
            fmd_fork_node_t *breakpoint;
            fmd_vector_t *incoming_sa_intervals;
            fmd_imt_query(fmd->r2r_tree, c_0_lo - 1, c_0_hi - 2, &incoming_sa_intervals);
            int_t no_forks_to_add = max_forks == -1 ? incoming_sa_intervals->size : MIN2(max_forks - forks->size, incoming_sa_intervals->size);
            if(no_forks_to_add) {
                breakpoint = fmd_fork_node_init(fork,
                                                fork->sa_lo, fork->sa_hi,
                                                fork->pos,
                                                BKPT);
                breakpoint->parent = (struct fmd_fork_node_t *)fork;
                fmd_fork_node_t *folded = fmd_fork_node_init(breakpoint,
                                                             fork->sa_lo, fork->sa_hi,
                                                             fork->pos,
                                                             MAIN);
                fork = folded;
                forks->data[i] = fork;
            }
            for (int_t j = 0; j < no_forks_to_add; j++) {
                fmd_imt_interval_t *interval = incoming_sa_intervals->data[j];
                fmd_fork_node_t *new_fork = fmd_fork_node_init(breakpoint,
                                                               V+1+interval->lo, V+2+interval->hi,
                                                               fork->pos,
                                                               FALT);
#pragma omp critical(forks_append)
                {
                    fmd_vector_append(new_forks, new_fork);
                }
            }
            fmd_vector_free(incoming_sa_intervals);
        }
    }
    /**********************************************************************
    * Step 2 - Merge phase: merge overlapping ranges on the vertices
    **********************************************************************/
    fmd_vector_t *merged;
    fmd_fmd_compact_forks(fmd, new_forks, &merged);
    fmd_vector_free(new_forks);
    /**********************************************************************
    * Step 3 - Advance Phase: advance each fork once
    **********************************************************************/
    fmd_vector_t *next_iter_forks;
    fmd_vector_init(&next_iter_forks, forks->size + merged->size, &fmd_fstruct_fork_node);
    // advance and filter previous queries
#pragma omp parallel for default(none) shared(forks, fmd, partial_matches, next_iter_forks, string, t)
    for(int_t i = 0; i < forks->size; i++) {
        fmd_fork_node_t *fork = forks->data[i];
        fmd_fmd_advance_fork(fmd, fork, string);
        if(fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fmd_fork_node_t *dead_node = fmd_fork_node_init(fork,
                                                            fork->sa_lo, fork->sa_hi,
                                                            fork->pos,
                                                            DEAD);
            fork = dead_node;
#pragma omp critical(partial_matches_append)
            {
                fmd_vector_append(*partial_matches, fork);
            }
        } else {
            if(*t==1) {
                fork->type = LEAF;
            }
#pragma omp critical(next_iter_queries_append)
            {
                fmd_vector_append(next_iter_forks, fork);
            }
        }
    }
    // advance and filter next forks
#pragma omp parallel for default(none) shared(merged,V,fmd,string,partial_matches,next_iter_forks, t)
    for (int_t i = 0; i < merged->size; i++) {
        fmd_fork_node_t *fork = merged->data[i];
        fmd_fmd_advance_fork(fmd, fork, string);
        if (fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fmd_fork_node_t *dead_node = fmd_fork_node_init(fork,
                                                            fork->sa_lo, fork->sa_hi,
                                                            fork->pos,
                                                            DEAD);
            fork = dead_node;
#pragma omp critical(partial_matches_append)
            {
                fmd_vector_append(*partial_matches, fork);
            }
        } else {
            fmd_fork_node_t *folded = fmd_fork_node_init(fork,
                                                         fork->sa_lo, fork->sa_hi,
                                                         fork->pos,
                                                         MAIN);
            if(*t==1) {
                folded->type = LEAF;
            }
            fork = folded;
#pragma omp critical(next_iter_queries_append)
            {
                fmd_vector_append(next_iter_forks, fork);
            }
        }
    }
    fmd_vector_free_disown(merged);
    fmd_vector_free_disown(forks);

    /************************************************************
    * Double compaction
    ************************************************************/
    /*
    fmd_vector_sort(next_iter_forks);
    fmd_vector_t *compacted;
    fmd_fmd_compact_forks(fmd, next_iter_forks, &compacted);
    fmd_vector_free(next_iter_forks);
    next_iter_forks = compacted;
    */
    /************************************************************
    * Beta feature, not yet benched!
    ************************************************************/
    *cur_forks = next_iter_forks;
    --(*t);
}

void fmd_fmd_query_find_bootstrapped(fmd_fmd_t *fmd, fmd_vector_t *bootstrap, int_t bootstrap_depth, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends) {
    fmd_vector_t *forks = bootstrap;
    fmd_vector_t *partial_matches;
    fmd_vector_init(&partial_matches, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);

    int_t t = bootstrap_depth; // stores the position last matched
    while(forks->size && t > 0) {
        fmd_fmd_query_find_step(fmd, string, max_forks, &t, &forks, &partial_matches);
    }
    *paths = forks;
    *dead_ends = partial_matches;
}

void fmd_fmd_query_find(fmd_fmd_t *fmd, fmd_fmd_cache_t *cache, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends) {
    fmd_vector_t *forks = NULL;
    int_t bootstrap_depth = -1;
    if(!cache) {
        // no cache: default bootstrap
        //data structures to keep track of forks
        fmd_vector_init(&forks, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);

        // create the initial fork
        int_t init_lo = 0;//1 + V * (2+fmd_ceil_log2(V));
        int_t init_hi = fmd->graph_fmi->no_chars;
        fmd_fork_node_t *root = fmd_fork_node_init(NULL,
                                                   init_lo, init_hi,
                                                   string->size - 1,
                                                   ROOT);
        // fold the root fork onto itself before advancing
        fmd_fork_node_t *main = fmd_fork_node_copy(root);
        main->parent = (struct fmd_fork_node_t *) root;
        main->type = MAIN;

        // now advance the root once before starting to fork
        fmd_fmd_advance_fork(fmd, main, string);

        if (main->sa_hi <= main->sa_lo) {
            fmd_vector_t *matches;
            fmd_vector_init(&matches, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
            *paths = matches;
            main->type = DEAD;
            fmd_vector_t *partial_matches;
            fmd_vector_init(&partial_matches, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
            fmd_vector_append(partial_matches, main);
            *dead_ends = partial_matches;
            fmd_vector_free(forks);
            return;
        }
        fmd_vector_append(forks, main);
        bootstrap_depth = string->size-1;
    } else {
        // cache lookup to get a bootstrap
        fmd_vector_t *bootstrap;
        void *lookup;
        int_t cached_suffix_start = string->size <= cache->depth ? 0 : string->size - cache->depth;
        fmd_string_t *cached_suffix;
        fmd_string_substring(string, cached_suffix_start, string->size, &cached_suffix);
        bootstrap_depth = cached_suffix_start;
        fmd_table_lookup(cache->tables[cached_suffix->size-1], cached_suffix, &lookup);
        fmd_vector_t *cached_forks = lookup;
        fmd_vector_init(&bootstrap, cached_forks->size, &fmd_fstruct_fork_node);
        #pragma omp parallel for default(none) shared(cached_forks, bootstrap, bootstrap_depth)
        for(int_t i = 0; i < cached_forks->size; i++) {
            fmd_fork_node_t *cached_fork = cached_forks->data[i];
            fmd_fork_node_type_t type = bootstrap_depth ? MAIN : LEAF;
            fmd_fork_node_t *fork = fmd_fork_node_init(cached_fork,
                                                       cached_fork->sa_lo,
                                                       cached_fork->sa_hi,
                                                       bootstrap_depth-1,
                                                       type);
            bootstrap->data[i] = fork;
        }
        bootstrap->size = cached_forks->size;
        forks = bootstrap;
    }
    fmd_fmd_query_find_bootstrapped(fmd, forks, bootstrap_depth, string, max_forks, paths, dead_ends);
}

void fmd_fmd_compact_forks(fmd_fmd_t *fmd, fmd_vector_t *forks, fmd_vector_t **merged_forks) {
    fmd_vector_sort(forks);
    fmd_vector_t *merged;
    fmd_vector_init(&merged, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    if(forks->size) {
        fmd_fork_node_t *cur_fork = fmd_fork_node_copy(forks->data[0]);
        fmd_vector_append(merged, cur_fork);
        for (int_t i = 1; i < forks->size; i++) {
            fmd_fork_node_t *next_fork = forks->data[i];
            if(cur_fork->sa_hi >= next_fork->sa_lo) {
                cur_fork->sa_hi = next_fork->sa_hi > cur_fork->sa_hi ? next_fork->sa_hi : cur_fork->sa_hi;
            } else {
                cur_fork = fmd_fork_node_copy(next_fork);
                fmd_vector_append(merged, cur_fork);
            }
        }
    }
    *merged_forks = merged;
}

void fmd_fmd_topologise_fork(fmd_fork_node_t *fork, fmd_string_t *query, fmd_match_chain_t **chain) {
    fmd_match_chain_t *match_chain;
    fmd_fmd_match_chain_init(&match_chain);
    fmd_fork_node_t *cur_fork = fork;
    int_t start = 0;
    int_t end = -1;
    int_t cur_sa_lo = -1;
    int_t cur_sa_hi = -1;
    bool terminate = false;
    while(!terminate) {
        switch (cur_fork->type) {
            case LEAF:{
                cur_sa_lo = cur_fork->sa_lo;
                cur_sa_hi = cur_fork->sa_hi;
                break;
            }
            case CACH:
            case ROOT: {
                end = query->size;
                fmd_string_t *subs;
                fmd_string_substring(query, start, end, &subs);
                fmd_fmd_match_chain_append(match_chain, subs, cur_sa_lo, cur_sa_hi);
                start = end;
                terminate = true;
                break;
            }
            case FALT: {
                end = ((fmd_fork_node_t*)(cur_fork->parent))->pos+1;
                fmd_string_t *subs;
                fmd_string_substring(query, start, end, &subs);
                fmd_fmd_match_chain_append(match_chain, subs, cur_sa_lo, cur_sa_hi);
                start = end;
                cur_fork = (fmd_fork_node_t*)cur_fork->parent;
                cur_sa_lo = cur_fork->sa_lo;
                cur_sa_hi = cur_fork->sa_hi;
                break;
            }
            default: {
                break;
            }
        }
        cur_fork = (fmd_fork_node_t*)cur_fork->parent;
        if(!cur_fork) terminate = true;
    }
    *chain = match_chain;
}

void fmd_fmd_topologise_forks(fmd_string_t *query, fmd_vector_t *input_matches, fmd_vector_t **match_chains, int_t *count) {
    if(!query || !input_matches->size) {
        fmd_vector_init(match_chains, 1, &fmd_fstruct_match_list);
    }
    fmd_vector_t *matches;
    fmd_vector_init(&matches, input_matches->size, &fmd_fstruct_match_list);

    int_t total_count = 0;

    for(int_t i = 0; i < input_matches->size; i++) {
        fmd_match_chain_t *match_chain;
        fmd_fork_node_t *fork = input_matches->data[i];
        total_count += fork->sa_hi - fork->sa_lo;
        fmd_fmd_topologise_fork(fork, query, &match_chain);
        fmd_vector_append(matches, match_chain);
    }
    *match_chains = matches;
    *count = total_count;
}

void fmd_fmd_topologise_forks_free(fmd_vector_t *match_lists) {
    if(!match_lists) return;
    fmd_vector_free(match_lists);
}

void fmd_decoded_match_init(fmd_decoded_match_t **dec, vid_t vid, int_t offset) {
    fmd_decoded_match_t *d = calloc(1, sizeof(fmd_decoded_match_t));
    d->vid = vid;
    d->offset = offset;
    *dec = d;
}

int fmd_decoded_match_comp(fmd_decoded_match_t *dec1, fmd_decoded_match_t *dec2) {
    if(!dec1 || !dec2) return -1;
    return dec1->vid == dec2->vid ? (int)(dec1->offset - dec2->offset) : (int)(dec1->vid - dec2->vid);
}

uint_t fmd_decoded_match_hash(fmd_decoded_match_t *dec) {
    return prm_hash_f((void*)dec->vid) ^ prm_hash_f((void*)dec->offset);
}

void fmd_decoded_match_free(fmd_decoded_match_t *dec) {
    if(dec) {
        free(dec);
    }
}

fmd_decoded_match_t* fmd_decoded_match_copy(fmd_decoded_match_t *dec) {
    if(!dec) return NULL;
    fmd_decoded_match_t *d = calloc(1, sizeof(fmd_decoded_match_t));
    d->vid = dec->vid;
    d->offset = dec->offset;
    return d;
}

void fmd_fmd_decoder_init(fmd_fmd_decoder_t **dec, fmd_fmd_t *fmd) {
    fmd_fmd_decoder_t *d = calloc(1, sizeof(fmd_fmd_decoder_t));
    if(!d) {
        *dec = NULL;
        return;
    }
    d->fmd = fmd;
    fmd_vector_init(&d->vertex_bases, fmd->permutation->size, &prm_fstruct);

    int_t V = fmd->permutation->size;
    fmd_fmi_qr_t qr;
    qr.lo = 1;
    qr.hi = V + 1;
    fmd_vector_t *bases_permuted = fmd_fmi_sa(fmd->graph_fmi, &qr);
    for(int_t i = 0; i < V; i++) { // code below actually sorts the array :)
        d->vertex_bases->data[(vid_t)fmd->bwt_to_vid->data[i]] = bases_permuted->data[i];
    }
    d->vertex_bases->size = V;

    fmd_vector_free(bases_permuted);
    *dec = d;
}

void fmd_fmd_decoder_free(fmd_fmd_decoder_t *dec) {
    if(!dec) return;
    fmd_vector_free(dec->vertex_bases);
    // don't free dec->fmd
    free(dec);
}

void fmd_fmd_decoder_decode_one(fmd_fmd_decoder_t *dec, int_t sa_lo, int_t sa_hi, int_t matches_to_decode, fmd_vector_t **matches) {
    int_t no_to_decode = matches_to_decode == -1 ? sa_hi - sa_lo : MIN2(matches_to_decode, sa_hi - sa_lo);
    fmd_vector_t *m;
    fmd_vector_init(&m, no_to_decode, &fmd_fstruct_decoded_match);

    // decode the suffix array in the interval
    fmd_fmi_qr_t whole_rec;
    whole_rec.lo = sa_lo;
    whole_rec.hi = sa_lo + no_to_decode;
    fmd_vector_t *T = fmd_fmi_sa(dec->fmd->graph_fmi, &whole_rec);

    #pragma omp parallel for default(none) shared(sa_lo, sa_hi, dec, m ,T, no_to_decode)
    for(int_t i = 0; i < no_to_decode; i++) {
        // find the closest preceding vid in text space via binary search
        int_t lo = 0;
        int_t hi = dec->vertex_bases->size-1;
        int_t sa_val = (int_t)T->data[i]; // cache in local stack
        vid_t closest = 0;
        while(lo <= hi) {
            int_t mid = lo + (hi-lo)/2;
            int_t midpoint = (int_t)dec->vertex_bases->data[mid];
            if(midpoint < sa_val) {
                closest = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        // decode into the match
        fmd_decoded_match_t *match;
        fmd_decoded_match_init(&match, closest, T->data[i] - dec->vertex_bases->data[closest] - 1);
        m->data[i] = match;
    }
    m->size = no_to_decode;
    fmd_vector_free(T);
    fmd_vector_sort(m);
    *matches = m;
}

void fmd_fmd_decoder_decode_ends(fmd_fmd_decoder_t *dec, fmd_vector_t *matches, int_t max_matches, fmd_vector_t **decoded) {
    fmd_vector_t *all_decoded;
    fmd_vector_init(&all_decoded, matches->size, &fmd_fstruct_vector);
    int_t total_decoded = 0;
    int_t i;
    for(i = 0; i < matches->size; i++) {
        fmd_fork_node_t *fork = matches->data[i];
        int_t fork_size = fork->sa_hi - fork->sa_lo;
        int_t no_to_decode_from_fork = max_matches == -1 || (total_decoded + fork_size <= max_matches) ? fork_size : MIN2(max_matches - total_decoded, fork_size);
        fmd_vector_t *decoded_fork;
        fmd_fmd_decoder_decode_one(dec, fork->sa_lo, fork->sa_hi, no_to_decode_from_fork, &decoded_fork);
        all_decoded->data[i] = decoded_fork;
        total_decoded += no_to_decode_from_fork;
        if(max_matches != -1 && total_decoded >= max_matches) {
            ++i;
            break;
        }
    }
    all_decoded->size = i;
    *decoded = all_decoded;
}

void fmd_fmd_query_find_result_free(fmd_vector_t *paths, fmd_vector_t *dead_ends) {
    // keys are pointers, values are actual pointers
    fmd_tree_t *visited;
    fmd_tree_init(&visited, &prm_fstruct, &fmd_fstruct_fork_node);

    for(int_t i = 0; i < paths->size; i++) {
        fmd_fork_node_t *cur = paths->data[i];
        while(cur && cur->type != CACH) {
            bool inserted = fmd_tree_insert(visited, cur, cur);
            cur = (fmd_fork_node_t*)cur->parent;
        }
    }
    for(int_t i = 0; i < dead_ends->size; i++) {
        fmd_fork_node_t *cur = dead_ends->data[i];
        while(cur && cur->type != CACH) {
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

/*
void fmd_fmd_cache_init_helper(fmd_fmd_cache_t *cache, fmd_fmd_t *fmd, fmd_string_t *str, fmd_vector_t *forks, fmd_vector_t *partial_matches) {
    // generate extensions of the string
    fmd_vector_t *alphabet = fmd->graph_fmi->alphabet;
    for(int_t i = 5; i < alphabet->size; i++) { // first five characters are RESERVED.
        fmd_string_t *extension = fmd_string_copy(str);
        fmd_string_insert(extension, 0, (char_t)alphabet->data[i]);
        fmd_vector_t *cur_forks = fmd_vector_copy(forks);
        fmd_vector_init(&partial_matches, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
        int_t t = 1;
        fmd_fmd_query_find_step(fmd, extension, -1, &t, &cur_forks, &partial_matches);
        for(int_t j = 0; j < cur_forks->size; j++) {
            fmd_fork_node_t *fork = cur_forks->data[j];
            fork->pos = 0;
            fork->type = CACH;
        }
        fmd_table_insert(&cache->cache[extension->size-1], extension, cur_forks);
        if(extension->size < cache->depth) {
            fmd_fmd_cache_init_helper(cache, fmd, extension, cur_forks, partial_matches);
        }
    }
}
 */

void fmd_fmd_cache_init_step(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **cur_forks, fmd_vector_t **partial_matches) {
    fmd_vector_t *forks= *cur_forks;
    int_t V = fmd->permutation->size;
    /**********************************************************************
    * Step 1 - Forking phase: Fork each query at its current position
    **********************************************************************/
    fmd_vector_t *new_forks;
    fmd_vector_init(&new_forks, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
#pragma omp parallel for default(none) shared(forks, fmd, max_forks, new_forks, V)
    for (int_t i = 0; i < forks->size; i++) {
        fmd_fork_node_t *fork = forks->data[i];

        //printf("orig:%s\n", string->seq+1);
        //printf("ext:%s\n", string->seq);
        //printf("fork: %lld,%lld\n", fork->sa_lo, fork->sa_hi);
        //printf("\n");

        int_t c_0_lo, c_0_hi;
        fmd_fmd_fork_precedence_range(fmd, fork, fmd->c_0, &c_0_lo, &c_0_hi);
        if(c_0_lo < c_0_hi) {
            fmd_fork_node_t *breakpoint;
            fmd_vector_t *incoming_sa_intervals;
            fmd_imt_query(fmd->r2r_tree, c_0_lo - 1, c_0_hi - 2, &incoming_sa_intervals);
            if(incoming_sa_intervals->size) {
                breakpoint = fmd_fork_node_init(fork,
                                                fork->sa_lo, fork->sa_hi,
                                                fork->pos,
                                                BKPT);
                breakpoint->parent = (struct fmd_fork_node_t *)fork;
                fmd_fork_node_t *folded = fmd_fork_node_init(breakpoint,
                                                             fork->sa_lo, fork->sa_hi,
                                                             fork->pos,
                                                             CACH);
                fork = folded;
                forks->data[i] = fork;
            }
            for (int_t j = 0; j < incoming_sa_intervals->size; j++) {
                fmd_imt_interval_t *interval = incoming_sa_intervals->data[j];
                fmd_fork_node_t *new_fork = fmd_fork_node_init(breakpoint,
                                                               V+1+interval->lo, V+2+interval->hi,
                                                               fork->pos,
                                                               FALT);
#pragma omp critical(forks_append)
                {
                    fmd_vector_append(new_forks, new_fork);
                }
            }
            fmd_vector_free(incoming_sa_intervals);
        }
    }
    /**********************************************************************
    * Step 2 - Merge phase: merge overlapping ranges on the vertices
    **********************************************************************/
    fmd_vector_t *merged;
    fmd_fmd_compact_forks(fmd, new_forks, &merged);
    fmd_vector_free(new_forks);
    /**********************************************************************
    * Step 3 - Advance Phase: advance each fork once
    **********************************************************************/
    fmd_vector_t *next_iter_forks;
    fmd_vector_init(&next_iter_forks, forks->size + merged->size, &fmd_fstruct_fork_node);
    // advance and filter previous queries
#pragma omp parallel for default(none) shared(forks, fmd, partial_matches, next_iter_forks, string, t)
    for(int_t i = 0; i < forks->size; i++) {
        fmd_fork_node_t *fork = forks->data[i];
        fmd_fmd_advance_fork(fmd, fork, string);
        if(fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fmd_fork_node_t *dead_node = fmd_fork_node_init(fork,
                                                            fork->sa_lo, fork->sa_hi,
                                                            -1,
                                                            DEAD);
            fork = dead_node;
#pragma omp critical(partial_matches_append)
            {
                fmd_vector_append(*partial_matches, fork);
            }
        } else {
            fork->type = CACH;
            fork->pos = 0; // rewind
#pragma omp critical(next_iter_queries_append)
            {
                fmd_vector_append(next_iter_forks, fork);
            }
        }
    }
    // advance and filter next forks
#pragma omp parallel for default(none) shared(merged,V,fmd,string,partial_matches,next_iter_forks, t)
    for (int_t i = 0; i < merged->size; i++) {
        fmd_fork_node_t *fork = merged->data[i];
        fmd_fmd_advance_fork(fmd, fork, string);
        if (fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fmd_fork_node_t *dead_node = fmd_fork_node_init(fork,
                                                            fork->sa_lo, fork->sa_hi,
                                                            -1,
                                                            DEAD);
            fork = dead_node;
#pragma omp critical(partial_matches_append)
            {
                fmd_vector_append(*partial_matches, fork);
            }
        } else {
            fmd_fork_node_t *folded = fmd_fork_node_init(fork,
                                                         fork->sa_lo, fork->sa_hi,
                                                         0, // rewind
                                                         CACH);
            fork = folded;
#pragma omp critical(next_iter_queries_append)
            {
                fmd_vector_append(next_iter_forks, fork);
            }
        }
    }
    fmd_vector_free_disown(merged);
    fmd_vector_free_disown(forks);
    *cur_forks = next_iter_forks;
}

void fmd_fmd_cache_init_helper_trav(void *key, void *value, void *params) {
    fmd_fmd_cache_t *cache = (fmd_fmd_cache_t*)params;
    fmd_string_t *base = (fmd_string_t*)key;
    fmd_vector_t *matching_forks = (fmd_vector_t*)value;

    fmd_vector_t *alphabet = cache->fmd->graph_fmi->alphabet;
    // for all extensions, advance the forks once
    // set functions of the original list to copy only the references
    matching_forks->f = &prm_fstruct;
    for(int i = FMD_FMD_NO_RESERVED_CHARS; i < alphabet->size; i++) {
        // get extension
        fmd_string_t *extension;
        char_t ch = (char_t)alphabet->data[i];
        fmd_string_init(&extension, base->size+1);
        fmd_string_append(extension, ch);
        fmd_string_concat_mut(extension, base);
        // get initial list of forks
        fmd_vector_t *ref_forks;
        fmd_vector_init(&ref_forks, matching_forks->size, &prm_fstruct);
        ref_forks->size = matching_forks->size;
        for(int_t j = 0; j < matching_forks->size; j++) {
            fmd_fork_node_t *fork = matching_forks->data[j];
            ref_forks->data[j] = fmd_fork_node_init(fork,
                                                    fork->sa_lo, fork->sa_hi,
                                                    fork->pos,
                                                    CACH);
        }
        // advance all forks once
        fmd_vector_t *partial_matches;
        fmd_vector_init(&partial_matches, FMD_VECTOR_INIT_SIZE, &prm_fstruct); // store only the references
        fmd_fmd_cache_init_step(cache->fmd, extension, &ref_forks, &partial_matches);
        // Tree pruning: (LITERALLY)
        // handle partial matches: construct a table of visited nodes from the forks to the root
        // then for each partial match, free until a parent is a visited node.
        fmd_tree_t *visited;
        fmd_tree_init(&visited, &prm_fstruct, &prm_fstruct);
        for(int_t j = 0; j < ref_forks->size; j++) {
            fmd_fork_node_t *cur = ref_forks->data[j];
            while(cur) {
                fmd_tree_insert(visited, cur, NULL);
                cur = (fmd_fork_node_t*)cur->parent;
            }
        }
        for(int_t j = 0; j < partial_matches->size; j++) {
            void* _;
            fmd_fork_node_t *cur = ref_forks->data[j];
            bool vis = fmd_tree_search(visited, cur, &_);
            while(cur && !vis) {
                fmd_fork_node_t *tmp = cur;
                cur = (fmd_fork_node_t*)cur->parent;
                fmd_fork_node_free(tmp);
                vis = fmd_tree_search(visited, cur, &_);
            }
        }
        fmd_tree_free(visited);
        fmd_vector_free_disown(partial_matches);
        // insert the extension into the cache, and the lists are pointer datatypes.
        ref_forks->f = &prm_fstruct;
        fmd_table_insert(cache->tables[extension->size-1], extension, ref_forks);
    }
}

void fmd_fmd_cache_init(fmd_fmd_cache_t **cache, fmd_fmd_t *fmd, int_t depth) {
    fmd_fmd_cache_t *c = calloc(1, sizeof(fmd_fmd_cache_t*));
    if(!c) {
        *cache = NULL;
        return;
    }
    c->tables = calloc(depth, sizeof(fmd_table_t*));
    c->depth = depth;
    c->fmd = fmd;
    // initialize the cache
    for(int_t i = 0; i < depth; i++) {
        fmd_table_init(&c->tables[i], FMD_HT_INIT_SIZE, &fmd_fstruct_string, &fmd_fstruct_vector);
    }
    // populate the cache
    int_t init_lo = 0;
    int_t init_hi = fmd->graph_fmi->no_chars;
    fmd_fork_node_t *root = fmd_fork_node_init(NULL,
                                               init_lo, init_hi,
                                               0,
                                               CACH);

    fmd_vector_t *alphabet = fmd->graph_fmi->alphabet;
    fmd_vector_t *partial_matches;
    fmd_vector_init(&partial_matches, FMD_VECTOR_INIT_SIZE, &fmd_fstruct_fork_node);
    for(int_t i = FMD_FMD_NO_RESERVED_CHARS; i < alphabet->size; i++) { // first five characters are RESERVED.
        char_t ch = (char_t)alphabet->data[i];
        fmd_string_t *str;
        fmd_string_init(&str, FMD_STRING_INIT_SIZE);
        fmd_string_insert(str, 0, ch);
        // fold the root fork
        fmd_fork_node_t *main = fmd_fork_node_copy(root);
        main->parent = (struct fmd_fork_node_t *) root;
        main->type = MAIN;
        main->pos = 0;
        fmd_fmd_advance_fork(fmd, main, str);
        main->pos = 0;
        fmd_vector_t *forks;
        fmd_vector_init(&forks, 1, &prm_fstruct);
        fmd_vector_append(forks, main);
        fmd_table_insert(c->tables[0], str, forks);
    }
    // extend each bucket
    c->fmd = fmd;
    for(int_t i = 0; i < depth-1; i++) {
        fmd_table_traverse(c->tables[i], c, fmd_fmd_cache_init_helper_trav);
    }
    *cache = c;
}

void fmd_fmd_cache_free(fmd_fmd_cache_t *cache) {
    return; // BIG TODO
}

bool fmd_fmd_advance_fork(fmd_fmd_t *fmd, fmd_fork_node_t *fork, fmd_string_t *pattern) {
    fmd_fmi_t *fmi = fmd->graph_fmi;
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    bool found = fmd_table_lookup(fmi->c2e, (void*)pattern->seq[fork->pos], &encoding);
    if(!found) {
        fork->sa_lo = 0;
        fork->sa_hi = 0;
        //fprintf(stderr,"[fmd_fmi_advance_query]: encoding not found in dictionary, query is NIL\n");
        return false;
    }
    count_t rank_lo_m_1 = fork->sa_lo ? fmd_fmi_rank(fmi,encoding, fork->sa_lo-1) : 0ull;
    count_t rank_hi_m_1 = fork->sa_hi ? fmd_fmi_rank(fmi,encoding, fork->sa_hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    fork->sa_lo = (int_t)(base + rank_lo_m_1);
    fork->sa_hi = (int_t)(base + rank_hi_m_1);
    --fork->pos;
    return true;
}

bool fmd_fmd_fork_precedence_range(fmd_fmd_t *fmd, fmd_fork_node_t *fork, char_t c, int_t *lo, int_t *hi) {
    fmd_fmi_t *fmi = fmd->graph_fmi;
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    bool found = fmd_table_lookup(fmi->c2e, (void*)c, &encoding);
    if(!found) {
        fork->sa_lo = 0;
        fork->sa_hi = 0;
        //fprintf(stderr,"[fmd_fmi_advance_query]: encoding not found in dictionary, query is NIL\n");
        return false;
    }
    count_t rank_lo_m_1 = fork->sa_lo ? fmd_fmi_rank(fmi,encoding, fork->sa_lo-1) : 0ull;
    count_t rank_hi_m_1 = fork->sa_hi ? fmd_fmi_rank(fmi,encoding, fork->sa_hi-1) : 0ull;
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
    fmd_bs_free(bs);
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