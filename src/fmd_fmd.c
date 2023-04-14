#include "fmd_fmd.h"

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
    fmd_table_t *cword_map;
    fmd_table_init(&cword_map, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    // if(!cword_map...)
    for(int_t i = 0; i < permutation->size; i++) {
        fmd_table_insert(cword_map, permutation->data[i], (void*)i);
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
        int_t pcode_idx = -1;
        fmd_table_lookup(cword_map, (void*)i, (void*)&pcode_idx);
        fmd_string_concat_mut(graph_encoding, (fmd_string_t*)(cwords->data[pcode_idx]));
    }
    // free irrelevant stuff
    fmd_vector_free(cwords);
    fmd_table_free(cword_map);
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
        assert(graph_encoding->seq[sa[i]] == c_0);
    }
    for(int_t i = V+1; i <= 2*V; i++) {
        fmd_vector_append(c_1_bucket, (void*)sa[i]);
        // DEBUG PURPOSES
        assert(graph_encoding->seq[sa[i]] == c_1);
    }
    // at this point, sa is no longer needed as slices are extraced
    free(sa); // frees gigs of memory :)
    /******************************************************
    * Step 3b - Compute rank translation tables
    ******************************************************/
    // get bwt ranks to text ranks mappings
    fmd_vector_t *c_0_bwt_to_text, *c_1_bwt_to_text, *c_1_text_to_bwt;
    fmd_vector_argsort(&c_0_bwt_to_text, c_0_bucket);
    fmd_vector_argsort(&c_1_bwt_to_text, c_1_bucket);

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
    // free intermediate structures
    fmd_vector_free(c_0_bucket);
    fmd_vector_free(c_1_bucket);
    fmd_vector_free(c_1_bwt_to_text);
    fmd_table_free(c_1_inversion_table);
    // DEBUG PURPOSES
    printf("C_0 BWT to text:\n");
    for(int_t i = 0; i < c_0_bwt_to_text->size; i++) {
        printf("%d:%d\n",i,c_0_bwt_to_text->data[i]);
    }
    printf("C_1 text to BWT:\n");
    for(int_t i = 0; i < c_1_bwt_to_text->size; i++) {
        printf("%d:%d\n",i,c_1_text_to_bwt->data[i]);
    }
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
            fmd_vector_append(neighbors_bwt_idx, (void*)c_1_text_to_bwt->data[(vid_t)neighbors->data[j]]);
        }
        fmd_vector_t *bwt_neighbor_intervals;
        fmd_vector_init(&bwt_neighbor_intervals, neighbors->size, &fmd_fstruct_imt_interval);
        if(neighbors_bwt_idx->size) {
            // DEBUG PURPOSES
            printf("bwt-range of %d: ", vid);
            for(int_t k = 0; k < neighbors_bwt_idx->size; k++) {
                printf("%d ", neighbors_bwt_idx->data[k]);
            }
            printf("\n");
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
            // free the singular list
            free(neighbors_bwt_idx);
        }
        fmd_vector_append(kv_pairs, bwt_neighbor_intervals);
    }
    f->bwt_to_vid = c_0_bwt_to_text;
    fmd_vector_free(c_1_text_to_bwt);
    /******************************************************
    * Step 3d - Now, actually construct the tree
    ******************************************************/
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
            // if(!ok... // this check is not needed if everything goes ok in coding time
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
            }
            bool okc = fmd_fmi_advance_query(fmd->graph_fmi, query);
            if(!okc || query->lo == query->hi) {
                break;
            }
        }
        count += query->hi - query->lo;
        fmd_fmi_qr_free(query);
    }
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
            // if(!ok... // this check is not needed if everything goes ok in coding time
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
    return locs;
}

void fmd_fmd_query_locate_paths(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **paths, fmd_vector_t **dead_ends) {
    // walk root count
    fmd_vector_t *leaves;
    fmd_vector_init(&leaves, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_vector_t *graveyard;
    fmd_vector_init(&graveyard, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
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
            // if(!ok... // this check is not needed if everything goes ok in coding time
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
    *paths = leaves;
    *dead_ends = graveyard;
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
    bool found = fmd_table_lookup(fmi->c2e, c, &encoding);
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
        //printf("%s\n", codeword);
        fmd_string_t *copy;
        fmd_string_init_cstr(&copy, codeword);
        fmd_vector_append(rval, copy);
    }
    return rval;
}