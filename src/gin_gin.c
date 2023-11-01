/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_gin.c is part of gin
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

#include "gin_gin.h"
#ifdef GIN_SDSL
#include "../wrapper/sdsl_wrapper.h"
#endif

gin_fork_node_t *gin_fork_node_init(int_t sa_lo, int_t sa_hi,
                                    int_t pos, gin_fork_node_type_t type) {
    gin_fork_node_t *fn = calloc(1, sizeof(gin_fork_node_t));
    fn->sa_lo = sa_lo;
    fn->sa_hi = sa_hi;
    fn->pos = pos;
    fn->type = type;
    return fn;
}

void gin_fork_node_free(gin_fork_node_t *node) {
    // frees only one node
    if(node) free(node);
}

gin_fork_node_t *gin_fork_node_copy(gin_fork_node_t *node) {
    if(!node) return NULL;
    gin_fork_node_t *copy = calloc(1, sizeof(gin_fork_node_t));
    if(!copy) {
        free(copy);
        return NULL;
    }
    memcpy(copy, node, sizeof(gin_fork_node_t));
    return copy;
}
uint_t gin_fork_node_hash(gin_fork_node_t *node) {
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

int gin_fork_node_comp_exact(gin_fork_node_t *n1, gin_fork_node_t *n2) {
    return memcmp(n1, n2, sizeof(gin_fork_node_t));
}
int gin_fork_node_comp(gin_fork_node_t *n1, gin_fork_node_t *n2) {
    return n1->sa_lo < n2->sa_lo ? -1 : (n1->sa_lo > n2->sa_lo);
}

void gin_gin_init(gin_gin_t** gin, gin_graph_t *graph, gin_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate) {
    gin_gin_t *f = calloc(1, sizeof(gin_gin_t));
    if (!f || !graph) {
        *gin = NULL;
        return;
    }
    f->c_0 = c_0;
    f->c_1 = c_1;
    /**************************************************************************
    * Step 0 - If no permutation is provided, generate an identity permutation
    **************************************************************************/
    bool permutation_alloc = false;
    if (!permutation) {
        permutation_alloc = true;
        gin_vector_init(&permutation, graph->vertex_list->size, &prm_fstruct);
        for (int_t i = 0; i < graph->vertex_list->size; i++) {
            gin_vector_append(permutation, (void *) i);
        }
    }
    assert(permutation->size == graph->vertex_list->size);
    assert(c_0 < c_1);
    f->permutation = gin_vector_copy(permutation);
    if (permutation_alloc) {
        gin_vector_free(permutation);
    }
    /**************************************************************************
    * Step 1 - Compute the graph encoding and the permutation encodings
    **************************************************************************/
    /******************************************************
    * Step 1a - Generate fixed length binary encodings
    ******************************************************/
    gin_vector_t *cwords = gin_gin_init_pcodes_fixed_binary_helper(GIN_GIN_DEFAULT_a_0, GIN_GIN_DEFAULT_a_1,
                                                                   graph->vertex_list->size);
    /******************************************************
    * Step 1b - Map the inverse of perm. to codewords
    ******************************************************/
    gin_vector_t *inverse_permutation;
    gin_vector_init(&inverse_permutation, f->permutation->size, &prm_fstruct);
    for (int_t i = 0; i < f->permutation->size; i++) {
        inverse_permutation->data[(int_t) f->permutation->data[i]] = (void *) i;
    }
    /******************************************************
    * Step 1c - Finally, create the graph encoding
    ******************************************************/
    // precompute size requirements
    int_t V = graph->vertex_list->size;
    int_t total_label_len = 0;
    for (int_t i = 0; i < V; i++) {
        total_label_len += ((gin_vertex_t *) graph->vertex_list->data[i])->label->size;
    }
    // concatenate everything
    gin_string_t *graph_encoding;
    gin_string_init(&graph_encoding, total_label_len + (2 + gin_ceil_log2(V)) * V);
    for (int_t i = 0; i < V; i++) {
        gin_string_append(graph_encoding, c_0);
        gin_string_concat_mut(graph_encoding, ((gin_vertex_t *) graph->vertex_list->data[i])->label);
        gin_string_append(graph_encoding, c_1);
        int_t pcode_idx = (int_t) inverse_permutation->data[i];
        gin_string_concat_mut(graph_encoding, (gin_string_t *) (cwords->data[pcode_idx]));
    }
    // free obsolete stuff
    /**************************************************************************
    * Step 2 - Compute the suffix array of the encoding and build an FMI
    **************************************************************************/
    gin_vector_free(cwords);
#ifdef GIN_SDSL
    f->graph_fmi = csa_wt_build(graph_encoding->seq, graph_encoding->size);
    csa_wt_populate_alphabet(f->graph_fmi, &f->alphabet, &f->alphabet_size);
    f->no_chars = csa_wt_bwt_length(f->graph_fmi);
#else
    gin_fmi_t *fmi;
    int64_t *sa = calloc(graph_encoding->size+1, sizeof(uint64_t));
    // if(!sa...)
    divsufsort64((sauchar_t*)graph_encoding->seq, (saidx64_t*)sa, graph_encoding->size+1);
    gin_fmi_init_with_sa(&fmi,graph_encoding,sa,rank_sample_rate,isa_sample_rate);
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
#endif
    /**************************************************************************
    * Step 3 - Compute the query range translation index
    **************************************************************************/
    // sa[0] -> suffix starting with \0, essentially the terminator
    // sa[1]   ... sa[V]  : range of c_0
    // sa[V+1] ... sa[2V] : range of c_1
    gin_vector_t *c_0_bucket, *c_1_bucket;
    gin_vector_init(&c_0_bucket, V, &prm_fstruct);
    gin_vector_init(&c_1_bucket, V, &prm_fstruct);
#ifdef GIN_SDSL
    csa_wt_sa(f->graph_fmi, (uint64_t *) c_0_bucket->data, 1, V);
    csa_wt_sa(f->graph_fmi, (uint64_t *) c_1_bucket->data, V + 1, 2 * V);
    c_0_bucket->size = V;
    c_1_bucket->size = V;
#else
    for(int_t i = 1; i <= V; i++) {
        //c_0_bucket->data[i] = (void*)sa[i];
        gin_vector_append(c_0_bucket, (void*)sa[i]);
        // DEBUG PURPOSES
        //assert(graph_encoding->seq[sa[i]] == c_0);
    }
    for(int_t i = V+1; i <= 2*V; i++) {
        gin_vector_append(c_1_bucket, (void*)sa[i]);
        // DEBUG PURPOSES
        //assert(graph_encoding->seq[sa[i]] == c_1);
    }
    // at this point, sa is no longer needed as slices are extracted

    free(sa);
#endif
    gin_string_free(graph_encoding);
    /******************************************************
    * Step 3b - Compute rank translation tables
    ******************************************************/
    // get bwt ranks to text ranks mappings
    gin_vector_t *c_0_bwt_to_text_tmp, *c_1_bwt_to_text, *c_1_text_to_bwt;
    gin_vector_argsort(&c_0_bwt_to_text_tmp, c_0_bucket);
    gin_vector_argsort(&c_1_bwt_to_text, c_1_bucket);

    // DEBUG PURPOSES
    gin_vector_t *c_0_bwt_to_text;
    gin_vector_argsort(&c_0_bwt_to_text, c_0_bwt_to_text_tmp);
    gin_vector_free(c_0_bwt_to_text_tmp);
    // END OF DEBUG PURPOSES

    // now, invert the mapping of c_1 to get a mapping from text to bwt for c_1 ranks
    gin_vector_init(&c_1_text_to_bwt, V, &prm_fstruct);
    gin_table_t *c_1_inversion_table;
    gin_table_init(&c_1_inversion_table, GIN_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    for(int_t i = 0; i < V; i++) {
        gin_table_insert(c_1_inversion_table, (void*)c_1_bwt_to_text->data[i], (void*)i);
    }
    for(int_t i = 0; i < V; i++) {
        int_t val = -1;
        bool found = gin_table_lookup(c_1_inversion_table, (void*)i, (void*)&val);
        assert(found);
        gin_vector_append(c_1_text_to_bwt, (void*)val);
    }

    // free intermediate structures
    gin_vector_free(c_0_bucket);
    gin_vector_free(c_1_bucket);
    gin_vector_free(c_1_bwt_to_text);
    gin_table_free(c_1_inversion_table);
    /******************************************************
    * Step 3c - (key,value)s for the interval merge tree
    ******************************************************/
    // keys : integers from 1 to V (vid+1)
    // values : list of intervals
    gin_vector_t *kv_pairs;
    gin_vector_init(&kv_pairs, V, &gin_fstruct_vector); // vector holds interval lists
    for(int_t i = 0; i < V; i++) {
        int_t vid = (int_t)c_0_bwt_to_text->data[i];
        gin_vector_t *neighbors;
        bool found = gin_table_lookup(graph->incoming_neighbors, (void*)vid, (void*)&neighbors);
        assert(found);
        gin_vector_t *neighbors_bwt_idx;
        gin_vector_init(&neighbors_bwt_idx, neighbors->size, &prm_fstruct);
        for(int_t j = 0; j < neighbors->size; j++) {
            gin_vector_append(neighbors_bwt_idx, (void*)inverse_permutation->data[(vid_t)neighbors->data[j]]);
        }
        gin_vector_t *bwt_neighbor_intervals;
        gin_vector_init(&bwt_neighbor_intervals, neighbors->size, &gin_fstruct_imt_interval);
        if(neighbors_bwt_idx->size) {
            // sort and compact neighbors_bwt_idx into intervals
            gin_vector_sort(neighbors_bwt_idx);
            // compaction start
            int_t lo = (int_t) neighbors_bwt_idx->data[0];
            int_t hi = lo;
            for (int_t j = 1; j < neighbors_bwt_idx->size; j++) {
                int_t current = (int_t) neighbors_bwt_idx->data[j];
                if (current == hi + 1) {
                    hi = current;
                } else {
                    gin_imt_interval_t *interval;
                    gin_imt_interval_init(&interval, lo, hi);
                    gin_vector_append(bwt_neighbor_intervals, interval);
                    lo = hi = current;
                }
            }
            // don't forget to add the last interval
            gin_imt_interval_t *interval;
            gin_imt_interval_init(&interval, lo, hi);
            gin_vector_append(bwt_neighbor_intervals, interval);
        }
        // free the singular list
        gin_vector_free(neighbors_bwt_idx);
        gin_vector_append(kv_pairs, bwt_neighbor_intervals);
    }
    f->bwt_to_vid = c_0_bwt_to_text;
    gin_vector_free(c_1_text_to_bwt);
    /******************************************************
    * Step 3d - Now, actually construct the tree
    ******************************************************/
    // END OF DEBUG PURPOSES
    gin_imt_t *inverval_merge_tree;
    gin_imt_init(&inverval_merge_tree, V, kv_pairs);
    f->r2r_tree = inverval_merge_tree;
    #ifdef GIN_ORACLE
    int_t *alphabet;
    int_t alphabet_size;
    int_t *vertex_last_char_enc;
    #ifdef GIN_SDSL
    alphabet = f->alphabet;
    alphabet_size = f->alphabet_size;
    vertex_last_char_enc = calloc(V, sizeof(int_t));
    csa_wt_bwt(f->graph_fmi, (uint64_t*)vertex_last_char_enc, V+1, 2*V);
    int_t *c2e = calloc(GIN_MAX_ALPHABET_SIZE, sizeof(int_t));
    int_t t = 0;
    for(int_t i = 0; i < alphabet_size; i++) {
        c2e[alphabet[i]] = t;
        t++;
    }
    for(int_t i = 0; i < V; i++) {
        vertex_last_char_enc[i] = c2e[vertex_last_char_enc[i]];
    }
    free(c2e);
    #else
    alphabet = f->graph_fmi->alphabet;
    alphabet_size = f->graph_fmi->alphabet_size;
    vertex_last_char_enc = calloc(V, sizeof(int_t));
    for(int_t i = 0; i < V; i++) {
        vertex_last_char_enc[i] = (int_t)gin_fmi_get(f->graph_fmi, V+1+i);
    }
    #endif
    gin_oimt_init(f->r2r_tree, vertex_last_char_enc, alphabet, alphabet_size, &f->oracle_r2r);
    free(vertex_last_char_enc);
    #endif
    // free the list structure storing interval lists, but not the lists themselves
    //kv_pairs->f = &prm_fstruct;
    gin_vector_free(kv_pairs);
    gin_vector_free(inverse_permutation);
    /******************************************************
    * Step 4 - Return
    ******************************************************/
    *gin = f;
}

void gin_gin_free(gin_gin_t *gin) {
    if(gin) {
        gin_vector_free(gin->permutation);
        gin_vector_free(gin->bwt_to_vid);
        free(gin->alphabet);
#ifdef GIN_SDSL
        csa_wt_free(gin->graph_fmi);
#else
        gin_fmi_free(gin->graph_fmi);
#endif
        gin_imt_free(gin->r2r_tree);
#if GIN_ORACLE
        gin_oimt_free(gin->oracle_r2r);
#endif
        free(gin);
    }
}

int gin_gin_comp(gin_gin_t *f1, gin_gin_t *f2) {
#ifdef GIN_SDSL
    int c1 = csa_wt_comp(f1->graph_fmi, f2->graph_fmi);
#else
    int c1 = gin_fmi_comp(f1->graph_fmi, f2->graph_fmi);
#endif
    int c2 = gin_imt_comp(f1->r2r_tree, f2->r2r_tree);
    int c3 = gin_vector_comp(f1->permutation, f2->permutation);
    int c4 = gin_vector_comp(f1->bwt_to_vid, f2->bwt_to_vid);
    int c5 = f1->c_0 == f2->c_0;
    int c6 = f1->c_1 == f2->c_1;
    return c1 == 0 && c2 == 0 && c3 == 0 && c4 == 0 && c5 && c6 ? 0 : -1;
}

#ifdef GIN_OMP
#include <omp.h>
#endif
void gin_gin_query_find_dfs(gin_gin_t *gin, gin_string_t *string, int_t max_forks, gin_vector_t **paths, gin_vector_t **dead_ends, int_t num_threads) {
    gin_vector_t *leaves;
    gin_vector_init(&leaves, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
    gin_vector_t *graveyard = NULL;
    gin_vector_init(&graveyard, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);

    // create initial task to fire
    int_t init_lo = 0;//1 + V * (2+gin_ceil_log2(V));
#ifdef GIN_SDSL
    int_t init_hi = gin->no_chars;
#else
    int_t init_hi = gin->graph_fmi->no_chars;
#endif

    gin_fork_node_t *root_fork = gin_fork_node_init(init_lo, init_hi,
                                                    string->size-1,
                                                    ROOT);
#ifdef GIN_OMP
    omp_set_num_threads((int)num_threads);
#endif
#pragma omp parallel default(none) shared(gin,root_fork,leaves,graveyard,max_forks,string)
    {
#pragma omp single nowait
        {
            gin_gin_query_find_dfs_process_fork(gin,root_fork,max_forks,string,leaves,graveyard);
        }
    }
#pragma omp taskwait
    *paths = leaves;
    *dead_ends = graveyard;
}

void gin_gin_query_find_dfs_process_fork(gin_gin_t *gin, gin_fork_node_t *fork, int_t max_forks, gin_string_t *pattern, gin_vector_t *exact_matches, gin_vector_t *partial_matches) {
    bool continue_flag = true;
#pragma omp critical(exact_match)
    {
        if(max_forks != -1 && exact_matches->size >= max_forks) {
            continue_flag = false;
        }
    }
    int_t V = gin->permutation->size;
    if (!continue_flag) {
        gin_fork_node_free(fork);
        return;
    }
    while (fork->pos > -1) {
        int_t c_0_lo, c_0_hi;
        bool okc = gin_gin_advance_fork(gin, fork, pattern);
        bool ok = gin_gin_fork_precedence_range(gin, fork, gin->c_0, &c_0_lo, &c_0_hi);
        if (fork->pos == -1) break;
        if (c_0_hi > c_0_lo) {
            // we have a walk having the current suffix of the query as a prefix
            gin_fork_node_t *royal_node = gin_fork_node_init(fork->sa_lo, fork->sa_hi,
                                                             fork->pos,
                                                             MAIN);
            gin_vector_t *incoming_sa_intervals;
            #ifdef GIN_ORACLE
            gin_oimt_query(gin->oracle_r2r, c_0_lo - 1, c_0_hi - 2, pattern->seq[fork->pos], max_forks, &incoming_sa_intervals);
            #else
            gin_imt_query(gin->r2r_tree, c_0_lo - 1, c_0_hi - 2, max_forks, &incoming_sa_intervals);
            #endif
            for (int_t i = 0; i < incoming_sa_intervals->size; i++) {
                gin_imt_interval_t *interval = incoming_sa_intervals->data[i];
                gin_fork_node_t *new_fork = gin_fork_node_init(V+1+interval->lo, V+2+interval->hi,
                                                               fork->pos, MAIN);
                // fire subtask
#pragma omp task default(none) shared(gin, new_fork, exact_matches, partial_matches, max_forks, pattern)
                gin_gin_query_find_dfs_process_fork(gin, new_fork, max_forks, pattern, exact_matches, partial_matches);
            }
            fork = royal_node;
            gin_vector_free(incoming_sa_intervals);
        }
        if (!okc || fork->sa_lo == fork->sa_hi) {
            fork->type = DEAD;
            break;
        }
    }
    if (fork->type != DEAD && fork->type != LEAF) {
        if (fork->sa_lo >= fork->sa_hi) {
            fork->type = DEAD;
#pragma omp critical(partial_match)
            {
                gin_vector_append(partial_matches, fork);
            }
        } else {
            fork->type = LEAF;
#pragma omp critical(exact_match)
            {
                if(max_forks == -1 || exact_matches->size < max_forks) {
                    gin_vector_append(exact_matches, fork);
                } else {
                    gin_vector_append(partial_matches, fork);
                }
            }
        }
    } else {
#pragma omp critical(partial_match)
        {
            gin_vector_append(partial_matches, fork);
        }
    }
    //gin_gin_qr_free(query);
}

void gin_gin_query_find_step(gin_gin_t *gin, gin_string_t *string, int_t max_forks, int_t *t, gin_vector_t **cur_forks, gin_vector_t **partial_matches, gin_gin_stats_t *stats) {
    gin_vector_t *forks= *cur_forks;
    int_t V = gin->permutation->size;
    /**********************************************************************
    * Step 1 - Forking phase: Fork each query at its current position
    **********************************************************************/
    gin_vector_t *new_forks;
    gin_vector_init(&new_forks, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
    #pragma omp parallel for default(none) shared(forks, gin, max_forks, new_forks, V)
    for (int_t i = 0; i < forks->size; i++) {
        gin_fork_node_t *fork = forks->data[i];
        int_t c_0_lo, c_0_hi;
        gin_gin_fork_precedence_range(gin, fork, gin->c_0, &c_0_lo, &c_0_hi);
        bool more_to_track = max_forks == -1 || max_forks > forks->size;
        if(more_to_track && c_0_lo < c_0_hi) {
            gin_vector_t *incoming_sa_intervals;
            #ifdef GIN_ORACLE
            gin_oimt_query(gin->oracle_r2r, c_0_lo - 1, c_0_hi - 2, string->seq[fork->pos], max_forks, &incoming_sa_intervals);
            #else
            gin_imt_query(gin->r2r_tree, c_0_lo - 1, c_0_hi - 2, max_forks, &incoming_sa_intervals);
            #endif
            int_t no_forks_to_add = max_forks == -1 ? incoming_sa_intervals->size : MIN2(max_forks - forks->size, incoming_sa_intervals->size);
            for (int_t j = 0; j < no_forks_to_add; j++) {
                gin_imt_interval_t *interval = incoming_sa_intervals->data[j];
                gin_fork_node_t *new_fork = gin_fork_node_init(V+1+interval->lo, V+2+interval->hi,
                                                               fork->pos,
                                                               MAIN);
                #pragma omp critical(forks_append)
                {
                    gin_vector_append(new_forks, new_fork);
                }
            }
            gin_vector_free(incoming_sa_intervals);
        }
    }
    /**********************************************************************
    * Step 2 - Merge phase: merge overlapping ranges on the vertices
    **********************************************************************/
    gin_vector_t *merged;
    gin_gin_compact_forks(gin, new_forks, &merged);
    gin_vector_free(new_forks);
    /**********************************************************************
    * Step 3 - Advance Phase: advance each fork once
    **********************************************************************/
    gin_vector_t *next_iter_forks;
    gin_vector_init(&next_iter_forks, forks->size + merged->size, &gin_fstruct_fork_node);
    // advance and filter previous queries
    #pragma omp parallel for default(none) shared(forks, gin, partial_matches, next_iter_forks, string, t)
    for(int_t i = 0; i < forks->size; i++) {
        gin_fork_node_t *fork = forks->data[i];
        gin_gin_advance_fork(gin, fork, string);
        if(fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fork->type = DEAD;
            #pragma omp critical(partial_matches_append)
            {
                gin_vector_append(*partial_matches, fork);
            }
        } else {
            if(*t==1) {
                fork->type = LEAF;
            }
            #pragma omp critical(next_iter_queries_append)
            {
                gin_vector_append(next_iter_forks, fork);
            }
        }
    }
    // advance and filter next forks
    #pragma omp parallel for default(none) shared(merged,V,gin,string,partial_matches,next_iter_forks, t)
    for (int_t i = 0; i < merged->size; i++) {
        gin_fork_node_t *fork = merged->data[i];
        gin_gin_advance_fork(gin, fork, string);
        if (fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fork->type = DEAD;
            #pragma omp critical(partial_matches_append)
            {
                gin_vector_append(*partial_matches, fork);
            }
        } else {
            #pragma omp critical(next_iter_queries_append)
            {
                gin_vector_append(next_iter_forks, fork);
            }
        }
    }
    stats->no_calls_to_advance_fork += forks->size + merged->size;
    stats->no_calls_to_precedence_range += forks->size;
    gin_vector_free_disown(merged);
    gin_vector_free_disown(forks);

    /************************************************************
    * Double compaction
    ************************************************************/
    /*
    //gin_vector_sort(next_iter_forks);
    gin_vector_t *compacted;
    gin_gin_compact_forks(gin, next_iter_forks, &compacted);
    gin_vector_free(next_iter_forks);
    next_iter_forks = compacted;
    */
    /************************************************************
    * Beta feature, not yet benched!
    ************************************************************/
    *cur_forks = next_iter_forks;
    --(*t);
}

void gin_gin_query_find_bootstrapped(gin_gin_t *gin, gin_vector_t *bootstrap, int_t bootstrap_depth, gin_string_t *string, int_t max_forks, gin_vector_t **paths, gin_vector_t **dead_ends, gin_gin_stats_t *stats) {
    gin_vector_t *forks = bootstrap;
    gin_vector_t *partial_matches;
    gin_vector_init(&partial_matches, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);

    int_t t = bootstrap_depth; // stores the position last matched
    while(forks->size && t > 0) {
        gin_gin_query_find_step(gin, string, max_forks, &t, &forks, &partial_matches, stats);
    }
    /* experimental */
    // compact forks one more time to prevent the reporting of duplicate matches
    //gin_vector_t *compacted;
    //gin_gin_compact_forks(gin, forks, &compacted);
    //gin_vector_free(forks);
    //forks = compacted;
    /* experimental */
    *paths = forks;
    *dead_ends = partial_matches;
    stats->no_matching_forks = forks->size;
    stats->no_partial_forks = partial_matches->size;
}

void gin_gin_query_find(gin_gin_t *gin,
                        gin_gin_cache_t *cache,
                        gin_string_t *string,
                        int_t max_forks,
                        gin_vector_t **paths,
                        gin_vector_t **dead_ends,
                        gin_gin_stats_t **stats) {
    gin_vector_t *forks = NULL;
    int_t bootstrap_depth = -1;
    *stats = calloc(1, sizeof(gin_gin_stats_t));
    if(!cache) {
        // no cache: default bootstrap
        //data structures to keep track of forks
        gin_vector_init(&forks, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);

        // create the initial fork
        int_t init_lo = 0;//1 + V * (2+gin_ceil_log2(V));
        int_t init_hi = gin->no_chars;
        gin_fork_node_t *root = gin_fork_node_init(init_lo, init_hi,
                                                   string->size - 1,
                                                   ROOT);
        // now advance the root once before starting to fork
        gin_gin_advance_fork(gin, root, string);

        if (root->sa_hi <= root->sa_lo) {
            gin_vector_t *matches;
            gin_vector_init(&matches, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
            *paths = matches;
            root->type = DEAD;
            gin_vector_t *partial_matches;
            gin_vector_init(&partial_matches, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
            gin_vector_append(partial_matches, root);
            *dead_ends = partial_matches;
            gin_vector_free(forks);
            ++(*stats)->no_calls_to_advance_fork;
            return;
        }
        gin_vector_append(forks, root);
        bootstrap_depth = string->size-1;
    } else {
        // cache lookup to get a bootstrap
        gin_vector_t *bootstrap;
        int_t cached_suffix_start = string->size <= cache->depth ? 0 : string->size - cache->depth;
        gin_string_t *cached_suffix;
        gin_string_substring(string, cached_suffix_start, string->size, &cached_suffix);
        bootstrap_depth = cached_suffix_start;
        gin_gin_cache_lookup(cache, cached_suffix, bootstrap_depth-1, max_forks, &bootstrap);
        gin_string_free(cached_suffix);
        forks = bootstrap;
    }
    gin_gin_query_find_bootstrapped(gin, forks, bootstrap_depth, string, max_forks, paths, dead_ends, *stats);
}

void gin_gin_compact_forks(gin_gin_t *gin, gin_vector_t *forks, gin_vector_t **merged_forks) {
    gin_vector_sort(forks);
    gin_vector_t *merged;
    gin_vector_init(&merged, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
    if(forks->size) {
        gin_fork_node_t *cur_fork = gin_fork_node_copy(forks->data[0]);
        gin_vector_append(merged, cur_fork);
        for (int_t i = 1; i < forks->size; i++) {
            gin_fork_node_t *next_fork = forks->data[i];
            if(cur_fork->sa_hi >= next_fork->sa_lo) {
                cur_fork->sa_hi = next_fork->sa_hi > cur_fork->sa_hi ? next_fork->sa_hi : cur_fork->sa_hi;
            } else {
                cur_fork = gin_fork_node_copy(next_fork);
                gin_vector_append(merged, cur_fork);
            }
        }
    }
    *merged_forks = merged;
}

void gin_decoded_match_init(gin_decoded_match_t **dec, vid_t vid, int_t offset) {
    gin_decoded_match_t *d = calloc(1, sizeof(gin_decoded_match_t));
    d->vid = vid;
    d->offset = offset;
    *dec = d;
}

int gin_decoded_match_comp(gin_decoded_match_t *dec1, gin_decoded_match_t *dec2) {
    if(!dec1 || !dec2) return -1;
    return dec1->vid == dec2->vid ? (int)(dec1->offset - dec2->offset) : (int)(dec1->vid - dec2->vid);
}

uint_t gin_decoded_match_hash(gin_decoded_match_t *dec) {
    return prm_hash_f((void*)dec->vid) ^ prm_hash_f((void*)dec->offset);
}

void gin_decoded_match_free(gin_decoded_match_t *dec) {
    if(dec) {
        free(dec);
    }
}

gin_decoded_match_t* gin_decoded_match_copy(gin_decoded_match_t *dec) {
    if(!dec) return NULL;
    gin_decoded_match_t *d = calloc(1, sizeof(gin_decoded_match_t));
    d->vid = dec->vid;
    d->offset = dec->offset;
    return d;
}

void gin_gin_decoder_init(gin_gin_decoder_t **dec, gin_gin_t *gin) {
    gin_gin_decoder_t *d = calloc(1, sizeof(gin_gin_decoder_t));
    if(!d) {
        *dec = NULL;
        return;
    }
    d->gin = gin;
    gin_vector_init(&d->vertex_bases, gin->permutation->size, &prm_fstruct);

    int_t V = gin->permutation->size;
    gin_fmi_qr_t qr;
    qr.lo = 1;
    qr.hi = V + 1;
#ifdef GIN_SDSL
    gin_vector_t *bases_permuted;
    gin_vector_init(&bases_permuted, V, &prm_fstruct);
    csa_wt_sa(gin->graph_fmi, (uint64_t*)bases_permuted->data, qr.lo, qr.hi-1);
    bases_permuted->size = qr.hi - qr.lo;
#else
    gin_vector_t *bases_permuted = gin_fmi_sa(gin->graph_fmi, &qr);
#endif
    for(int_t i = 0; i < V; i++) { // code below actually sorts the array :)
        d->vertex_bases->data[(vid_t)gin->bwt_to_vid->data[i]] = bases_permuted->data[i];
    }
    d->vertex_bases->size = V;

    gin_vector_free(bases_permuted);
    *dec = d;
}

void gin_gin_decoder_free(gin_gin_decoder_t *dec) {
    if(!dec) return;
    gin_vector_free(dec->vertex_bases);
    // don't free dec->gin
    free(dec);
}

void gin_gin_decoder_decode_one(gin_gin_decoder_t *dec, int_t sa_lo, int_t sa_hi, int_t matches_to_decode, gin_vector_t **matches) {
    int_t no_to_decode = matches_to_decode == -1 ? sa_hi - sa_lo : MIN2(matches_to_decode, sa_hi - sa_lo);
    gin_vector_t *m;
    gin_vector_init(&m, no_to_decode, &gin_fstruct_decoded_match);

    // decode the suffix array in the interval
    gin_fmi_qr_t whole_rec;
    whole_rec.lo = sa_lo;
    whole_rec.hi = sa_lo + no_to_decode;
    gin_vector_t *T;
#ifdef GIN_SDSL
    gin_vector_init(&T, no_to_decode, &prm_fstruct);
    csa_wt_sa(dec->gin->graph_fmi, (uint64_t*)T->data, sa_lo, sa_hi-1);
#else
    T = gin_fmi_sa(dec->gin->graph_fmi, &whole_rec);
#endif

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
        gin_decoded_match_t *match;
        gin_decoded_match_init(&match, closest, T->data[i] - dec->vertex_bases->data[closest] - 1);
        m->data[i] = match;
    }
    m->size = no_to_decode;
    gin_vector_free(T);
    gin_vector_sort(m);
    *matches = m;
}

void gin_gin_decoder_decode_ends(gin_gin_decoder_t *dec, gin_vector_t *matches, int_t max_matches, gin_vector_t **decoded) {
    gin_vector_t *all_decoded;
    gin_vector_init(&all_decoded, matches->size, &gin_fstruct_vector);
    int_t total_decoded = 0;
    int_t i;
    for(i = 0; i < matches->size; i++) {
        gin_fork_node_t *fork = matches->data[i];
        int_t fork_size = fork->sa_hi - fork->sa_lo;
        int_t no_to_decode_from_fork = max_matches == -1 || (total_decoded + fork_size <= max_matches) ? fork_size : MIN2(max_matches - total_decoded, fork_size);
        gin_vector_t *decoded_fork;
        gin_gin_decoder_decode_one(dec, fork->sa_lo, fork->sa_hi, no_to_decode_from_fork, &decoded_fork);
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

void gin_gin_cache_init_step(gin_gin_t *gin, gin_string_t *string, gin_vector_t **cur_forks, gin_vector_t **partial_matches) {
    gin_vector_t *forks= *cur_forks;
    int_t V = gin->permutation->size;
    /**********************************************************************
    * Step 1 - Forking phase: Fork each query at its current position
    **********************************************************************/
    gin_vector_t *new_forks;
    gin_vector_init(&new_forks, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
    for (int_t i = 0; i < forks->size; i++) {
        gin_fork_node_t *fork = forks->data[i];
        int_t c_0_lo, c_0_hi;
        gin_gin_fork_precedence_range(gin, fork, gin->c_0, &c_0_lo, &c_0_hi);
        if(c_0_lo < c_0_hi) {
            gin_vector_t *incoming_sa_intervals;
            gin_imt_query(gin->r2r_tree, c_0_lo - 1, c_0_hi - 2, -1, &incoming_sa_intervals);
            for (int_t j = 0; j < incoming_sa_intervals->size; j++) {
                gin_imt_interval_t *interval = incoming_sa_intervals->data[j];
                gin_fork_node_t *new_fork = gin_fork_node_init(V+1+interval->lo, V+2+interval->hi,
                                                               fork->pos, MAIN);
                gin_vector_append(new_forks, new_fork);
            }
            gin_vector_free(incoming_sa_intervals);
        }
    }
    /**********************************************************************
    * Step 2 - Merge phase: merge overlapping ranges on the vertices
    **********************************************************************/
    gin_vector_t *merged;
    gin_gin_compact_forks(gin, new_forks, &merged);
    gin_vector_free(new_forks);
    /**********************************************************************
    * Step 3 - Advance Phase: advance each fork once
    **********************************************************************/
    gin_vector_t *next_iter_forks;
    gin_vector_init(&next_iter_forks, forks->size + merged->size, &gin_fstruct_fork_node);
    // advance and filter previous queries
    for(int_t i = 0; i < forks->size; i++) {
        gin_fork_node_t *fork = forks->data[i];
        gin_gin_advance_fork(gin, fork, string);
        if(fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fork->type = DEAD;
            gin_vector_append(*partial_matches, fork);
        } else {
            fork->type = CACH;
            fork->pos = 0; // rewind
            gin_vector_append(next_iter_forks, fork);
        }
    }
    // advance and filter next forks
    for (int_t i = 0; i < merged->size; i++) {
        gin_fork_node_t *fork = merged->data[i];
        gin_gin_advance_fork(gin, fork, string);
        if (fork->sa_lo >= fork->sa_hi) { // query died while advancing
            fork->type = DEAD;
            gin_vector_append(*partial_matches, fork);
        } else {
            fork->type = CACH;
            fork->pos = 0; // rewind
            gin_vector_append(next_iter_forks, fork);
        }
    }
    gin_vector_free_disown(merged);
    gin_vector_free_disown(forks);
    *cur_forks = next_iter_forks;
}

void gin_gin_cache_init_helper_trav1(void *key, void *value, void *params) {
    gin_gin_cache_helper_p_t *p = params;
    gin_gin_cache_t *cache = p->cache;
    gin_gin_t *gin = p->gin;
    gin_table_t **cache_tables = p->cache_tables;
    gin_string_t *base = (gin_string_t*)key;
    gin_vector_t *matching_forks = (gin_vector_t*)value;
#ifdef GIN_SDSL
    int_t *alphabet = gin->alphabet;
    int_t alphabet_size = gin->alphabet_size;
#else
    int_t *alphabet = gin->graph_fmi->alphabet;
    int_t alphabet_size = gin->graph_fmi->alphabet_size;
#endif
    // for all extensions, advance the forks once
    // set functions of the original list to copy only the references
    for(int i = GIN_GIN_NO_RESERVED_CHARS; i < alphabet_size; i++) {
        // get extension
        gin_string_t *extension;
        char_t ch = (char_t)alphabet[i];
        gin_string_init(&extension, base->size+1);
        gin_string_append(extension, ch);
        gin_string_concat_mut(extension, base);
        // get initial list of forks
        gin_vector_t *ref_forks;
        gin_vector_init(&ref_forks, matching_forks->size, &gin_fstruct_fork_node);
        ref_forks->size = matching_forks->size;
        for(int_t j = 0; j < matching_forks->size; j++) {
            gin_fork_node_t *fork = matching_forks->data[j];
            ref_forks->data[j] = gin_fork_node_init(fork->sa_lo, fork->sa_hi,
                                                    fork->pos,
                                                    CACH);
        }
        // advance all forks once
        gin_vector_t *partial_matches;
        gin_vector_init(&partial_matches, GIN_VECTOR_INIT_SIZE, &gin_fstruct_fork_node);
        gin_gin_cache_init_step(gin, extension, &ref_forks, &partial_matches);
        gin_vector_free(partial_matches); // not necessary
        // insert the extension and its forks into the cache
        if(ref_forks->size) {
            /* experimental */
            // compact forks before inserting them into the table
            gin_vector_t *compacted;
            gin_gin_compact_forks(gin, ref_forks, &compacted);
            gin_vector_free(ref_forks);
            ref_forks = compacted;
            /* experimental */
            gin_table_insert(cache_tables[extension->size - 1], extension, ref_forks);
            ++cache->no_entries;
        } else {
            gin_string_free(extension);
            gin_vector_free(ref_forks);
        }
    }
}

void gin_gin_cache_init_helper_trav2(void* key, void* value, void* params) { //(*ftrav_kv)(void *key, void *value, void *p);
    gin_gin_cache_encode_p_t *p = params;
    gin_string_t *encoding = p->key_encoding;
    gin_vector_t *values = p->values;
    gin_string_append(encoding, GIN_GIN_DEFAULT_c_0);
    gin_string_concat_mut(encoding, (gin_string_t*)key);
    gin_vector_append(values, (gin_vector_t*)value);
}

void gin_gin_cache_init(gin_gin_cache_t **cache, gin_gin_t *gin, int_t depth) {
    /**********************************************************************
    * Step 0 - Initialize the main and helper data structures
    **********************************************************************/
    gin_gin_cache_t *c = calloc(1, sizeof(gin_gin_cache_t));
    if(!c) {
        *cache = NULL;
        return;
    }
    gin_table_t **cache_tables = calloc(depth, sizeof(gin_table_t*));
    c->depth = depth;
    for(int_t i = 0; i < depth; i++) {
        gin_table_init(&cache_tables[i], GIN_HT_INIT_SIZE, &gin_fstruct_string, &gin_fstruct_vector);
    }
    /**********************************************************************
    * Step 1 - Seed the cache with length 1 queries
    **********************************************************************/
    int_t init_lo = 0;
#ifdef GIN_SDSL
    int_t init_hi = gin->no_chars;
    int_t *alphabet = gin->alphabet;
    int_t alphabet_size = gin->alphabet_size;
#else
    int_t init_hi = gin->graph_fmi->no_chars;
    int_t *alphabet = gin->graph_fmi->alphabet;
    int_t alphabet_size = gin->graph_fmi->alphabet_size;
#endif
    for(int_t i = GIN_GIN_NO_RESERVED_CHARS; i < alphabet_size; i++) { // first five characters are RESERVED.
        char_t ch = (char_t)alphabet[i];
        gin_string_t *str;
        gin_string_init(&str, 1);
        gin_string_append(str, ch);
        gin_fork_node_t *root = gin_fork_node_init(init_lo, init_hi,
                                                   0,
                                                   MAIN);
        gin_gin_advance_fork(gin, root, str);
        if(root->sa_hi > root->sa_lo) { // only append to the table if the character exists
            root->pos = 0;
            gin_vector_t *forks;
            gin_vector_init(&forks, 1, &gin_fstruct_fork_node);
            gin_vector_append(forks, root);
            gin_table_insert(cache_tables[0], str, forks);
            ++c->no_entries;
        } else {
            gin_string_free(str);
            gin_fork_node_free(root);
        }
    }
    /**********************************************************************
    * Step 2 - Extend each seed recursively to matches until depth k
    **********************************************************************/
    gin_gin_cache_helper_p_t helper_params;
    helper_params.cache = c;
    helper_params.cache_tables = cache_tables;
    helper_params.gin = gin;
    // populate the tables of the cache
    for(int_t i = 0; i < depth-1; i++) {
        gin_table_traverse(cache_tables[i], &helper_params, gin_gin_cache_init_helper_trav1);
    }
    /**********************************************************************
    * Step 3 - Collapse table entries into a single string encoding with
    * the special character c_0
    **********************************************************************/
    gin_gin_cache_encode_p_t encode_params;
    encode_params.cache_tables = cache_tables;
    gin_string_init(&encode_params.key_encoding, GIN_STRING_INIT_SIZE);
    gin_vector_init(&encode_params.values, c->no_entries, &prm_fstruct); // stores pointers
    for(int_t i = 0; i < depth; i++) {
        gin_table_traverse(cache_tables[i], &encode_params, gin_gin_cache_init_helper_trav2);
    }
    /**********************************************************************
    * Step 4 - Construct the FM-Index over the special encoding and
    * rearrange the value list to match the ranks of the special character
    * c_0
    **********************************************************************/
#ifdef GIN_SDSL
    sdsl_csa *key_fmi = csa_wt_build(encode_params.key_encoding->seq, encode_params.key_encoding->size);
    uint64_t *sa = calloc(c->no_entries, sizeof(uint64_t));//printf("%s\n", encode_params.key_encoding->seq);
    csa_wt_sa(key_fmi, sa, 1, c->no_entries);
#else
    // first construct the suffix array
    int64_t *sa = calloc(encode_params.key_encoding->size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)encode_params.key_encoding->seq, (saidx64_t*)sa, encode_params.key_encoding->size+1);
    // then construct the FM-index
    gin_fmi_t *key_fmi;
    gin_fmi_init_with_sa(&key_fmi,
                         encode_params.key_encoding,
                         sa,
                         GIN_GIN_CACHE_FMI_DEFAULT_RANK_RATE,
                         encode_params.key_encoding->size+1);
#endif
    gin_string_free(encode_params.key_encoding);
    // now rearrange the items
    gin_vector_t *rearranged_values, *args;
    gin_vector_init(&rearranged_values, c->no_entries, &prm_fstruct);
#ifdef GIN_SDSL
    for(int_t i = 0; i < c->no_entries; i++) {
        rearranged_values->data[i] = (void*)sa[i];
    }
#else
    for(int_t i = 1; i <= c->no_entries; i++) {
        rearranged_values->data[i-1] = (void*)sa[i];
    }
#endif
    rearranged_values->size = c->no_entries;
    free(sa);
    gin_vector_argsort(&args, rearranged_values);
    // invert the argsort
    int_t no_value_bytes = 0;
    for(int_t i = 0; i < c->no_entries; i++) {
        gin_vector_t *interval_list = encode_params.values->data[i];
        rearranged_values->data[(int_t)args->data[i]] = interval_list;

        no_value_bytes += (GIN_GIN_CACHE_FORK_CARDINALITY_BIT_LENGTH >> 3) +
                          interval_list->size * (GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH >> 2);
    }
    gin_vector_free(args);
    /**********************************************************************
    * Step 4 - Write the value list and calculate byte offsets into the
    * buffer to read from
    **********************************************************************/
    gin_bs_t* value_bits;
    gin_bs_init_reserve(&value_bits, no_value_bytes / sizeof(word_t));
    word_t *word_offsets = calloc(c->no_entries, sizeof(word_t));
    int_t pos = 0;
    for(int_t i = 0; i < c->no_entries; i++) {
        word_offsets[i] = pos >> WORD_LOG_BITS;
        gin_vector_t *interval_list = rearranged_values->data[i];
        gin_bs_write_word(value_bits,
                          pos,
                          interval_list->size,
                          GIN_GIN_CACHE_FORK_CARDINALITY_BIT_LENGTH);
        pos += GIN_GIN_CACHE_FORK_CARDINALITY_BIT_LENGTH;
        for(int_t j = 0; j < interval_list->size; j++) {
            gin_fork_node_t *fork = interval_list->data[j];
            gin_bs_write_word(value_bits,
                              pos,
                              fork->sa_lo,
                              GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH);
            pos += GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH;
            gin_bs_write_word(value_bits,
                              pos,
                              fork->sa_hi,
                              GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH);
            pos += GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH;
        }
    }
    /**********************************************************************
    * Step 5 - Set values and clean up the data structures used
    **********************************************************************/
    gin_vector_free_disown(rearranged_values);
    for(int_t i = 0; i < depth; i++) {
        gin_table_free(cache_tables[i]);
    }
    free(cache_tables);
    uint_t no_words_value_bits;
    gin_bs_fit(value_bits, pos);
    gin_bs_detach(value_bits, &c->items, &no_words_value_bits);
    gin_bs_free(value_bits);
    gin_vector_free(encode_params.values);
    c->key_fmi = key_fmi;
    c->item_offsets = word_offsets;
    c->value_buffer_size_in_bits = (int_t)no_words_value_bits << WORD_LOG_BITS;
#ifdef GIN_SDSL
    c->key_fmi_size_in_bits = csa_wt_size_in_bytes(key_fmi) << 3;
#else
    c->key_fmi_size_in_bits = key_fmi->no_bits;
#endif
    *cache = c;
}

bool gin_gin_cache_advance_query(gin_gin_cache_t *cache, gin_gin_cache_qr_t *qr) {
#ifdef GIN_SDSL
    char c = qr->pattern->seq[qr->pos];
    count_t rank_lo = qr->lo ? csa_wt_rank(cache->key_fmi, qr->lo, c) : 0ull;
    count_t rank_hi = csa_wt_rank(cache->key_fmi, qr->hi, c);
    uint64_t base = csa_wt_char_sa_base(cache->key_fmi, c);
    qr->lo = (int_t)(base + rank_lo);
    qr->hi = (int_t)(base + rank_hi);
    --qr->pos;
    return qr->hi > qr->lo;
#else
    word_t encoding = cache->key_fmi->c2e[qr->pattern->seq[qr->pos]];
    count_t rank_lo_m_1 = qr->lo ? gin_fmi_rank(cache->key_fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? gin_fmi_rank(cache->key_fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = cache->key_fmi->char_counts[encoding];
    qr->lo = (int_t)(base + rank_lo_m_1);
    qr->hi = (int_t)(base + rank_hi_m_1);
    --qr->pos;
    return qr->hi > qr->lo;
#endif
}

bool gin_gin_cache_query_precedence_range(gin_gin_cache_t *cache, gin_gin_cache_qr_t *qr, char_t c, int_t *lo, int_t *hi) {
#ifdef GIN_SDSL
    count_t rank_lo = qr->lo ? csa_wt_rank(cache->key_fmi, qr->lo, c) : 0ull;
    count_t rank_hi = csa_wt_rank(cache->key_fmi, qr->hi, c);
    uint64_t base = csa_wt_char_sa_base(cache->key_fmi, c);
    *lo = (int_t)(base + rank_lo);
    *hi = (int_t)(base + rank_hi);
    return *hi > *lo;
#else
    word_t encoding = cache->key_fmi->c2e[c];
    count_t rank_lo_m_1 = qr->lo ? gin_fmi_rank(cache->key_fmi, encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? gin_fmi_rank(cache->key_fmi, encoding, qr->hi-1) : 0ull;
    uint64_t base = cache->key_fmi->char_counts[encoding];
    *lo = (int_t)(base + rank_lo_m_1);
    *hi = (int_t)(base + rank_hi_m_1);
    return *hi > *lo;
#endif
}

void gin_gin_cache_lookup(gin_gin_cache_t *cache, gin_string_t *string, int_t start_pos, int_t max_forks, gin_vector_t **cached_forks) {
    gin_vector_t *forks;
    gin_gin_cache_qr_t query;
    query.lo = 0;
#ifdef GIN_SDSL
    query.hi = csa_wt_bwt_length(cache->key_fmi);
#else
    query.hi = cache->key_fmi->no_chars+1;
#endif
    query.pos = string->size-1;
    query.pattern = string;
    while(query.pos > -1) {
        if(!gin_gin_cache_advance_query(cache, &query)) {
            gin_vector_init(&forks, 0, &gin_fstruct_fork_node);
            *cached_forks = forks;
            return;
        }
    }
    int_t lo, hi;
    if(!gin_gin_cache_query_precedence_range(cache, &query, GIN_GIN_DEFAULT_c_0, &lo, &hi)) {
        gin_vector_init(&forks, 0, &gin_fstruct_fork_node);
        *cached_forks = forks;
        return;
    }
    word_t start_offset = cache->item_offsets[lo-1];
    word_t *list = cache->items + start_offset + 1;
    word_t list_size = (int_t)cache->items[start_offset];
    word_t no_forks = max_forks == -1 ? list_size : (max_forks < list_size ? max_forks : list_size);
    gin_vector_init(&forks, (int_t)no_forks, &gin_fstruct_fork_node);
    forks->size = (int_t)no_forks;
    for(int_t i = 0; i < no_forks; i++) {
        int_t read_pos = i << 1;
        gin_fork_node_t *fork = gin_fork_node_init((int_t)list[read_pos],
                                                   (int_t)list[read_pos+1],
                                                   start_pos,
                                                   (start_pos == -1 ? LEAF : MAIN));
        forks->data[i] = fork;
    }
    *cached_forks = forks;
}

int_t gin_gin_cache_size(gin_gin_cache_t *cache) {
    if(!cache) return 0;
    int_t size = 0;
    size += GIN_GIN_CACHE_DEPTH_BIT_LENGTH;
    size += GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH;
    size += GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH;
    size += GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH;
    size += WORD_NUM_BITS * cache->no_entries;
    size += cache->key_fmi_size_in_bits;
    size += cache->value_buffer_size_in_bits;
    return (1 + ((size-1) >> 3));
}

void gin_gin_cache_serialize_to_buffer(gin_gin_cache_t *cache, unsigned char **buf_ret, uint64_t *buf_size_ret) {
    gin_bs_t *bs;
    gin_bs_init(&bs);
    uint_t widx = 0;
    /**************************************************************************
    * Step 0 - Write the header
    **************************************************************************/
    gin_bs_write_word(bs, widx, cache->depth, GIN_GIN_CACHE_DEPTH_BIT_LENGTH);
    widx += GIN_GIN_CACHE_DEPTH_BIT_LENGTH;
    gin_bs_write_word(bs, widx, cache->no_entries, GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH);
    widx += GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH;
    gin_bs_write_word(bs, widx, cache->key_fmi_size_in_bits, GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH);
    widx += GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH;
    gin_bs_write_word(bs, widx, cache->value_buffer_size_in_bits, GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH);
    widx += GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH;
    // word align
    widx = (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 1 - Write item offsets
    **************************************************************************/
    for (int_t i = 0; i < cache->no_entries; i++) {
        gin_bs_write_word(bs, widx, cache->item_offsets[i], WORD_NUM_BITS);
        widx += WORD_NUM_BITS;
    }
    // word align
    widx = (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 2 - Write the FMI bitstream
    **************************************************************************/

    uint8_t *fmi_buf = NULL;
    uint64_t fmi_buf_size = 0;
#ifdef GIN_SDSL
    csa_wt_to_buffer(cache->key_fmi, &fmi_buf, &fmi_buf_size);
#else
    gin_fmi_serialize_to_buffer(cache->key_fmi, &fmi_buf, &fmi_buf_size);
#endif
    int_t fmi_buf_size_in_bits = (int_t)fmi_buf_size << 3;
    gin_bs_write_word(bs, widx, (word_t)fmi_buf_size_in_bits, GIN_GIN_FMI_NO_BITS_BIT_LENGTH);
    widx += GIN_GIN_FMI_NO_BITS_BIT_LENGTH;
    int_t padded_no_bits = (1 + ((fmi_buf_size_in_bits-1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    int_t fmi_buf_no_padded_words = padded_no_bits >> WORD_LOG_BITS;
    word_t *fmi_buf_word = (word_t*)fmi_buf;
    for(int_t i = 0; i < fmi_buf_no_padded_words; i++) {
        gin_bs_write_word(bs, widx, fmi_buf_word[i], WORD_NUM_BITS);
        widx += WORD_NUM_BITS;
    }
    free(fmi_buf);
    // word align
    widx  = (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 3 - Write the item buffer
    **************************************************************************/
    int_t values_no_words = (int_t)(1 + ((cache->value_buffer_size_in_bits-1)>>WORD_LOG_BITS));
    for(int_t i = 0; i < values_no_words; i++) {
        gin_bs_write_word(bs, widx, (word_t)cache->items[i], WORD_NUM_BITS);
        widx += WORD_NUM_BITS;
    }
    /**************************************************************************
    * Step 4 - Return the buffer and free the other stuff
    **************************************************************************/
    word_t* buf;
    uint64_t no_words;
    gin_bs_fit(bs, widx);
    gin_bs_detach(bs, &buf, &no_words);
    gin_bs_free(bs);
    *buf_ret = (unsigned char*)buf;
    *buf_size_ret = no_words * sizeof(word_t);
}

void gin_gin_cache_serialize_from_buffer(gin_gin_cache_t **cachew, unsigned char *buf, uint64_t buf_size) {
    gin_bs_t *bs;
    uint_t ridx = 0;
    gin_gin_cache_t *cache;
    gin_bs_init_from_buffer(buf, (size_t)buf_size, &bs);
    cache = calloc(1, sizeof(gin_gin_cache_t));
    /**************************************************************************
    * Step 0 - Parse the header
    **************************************************************************/
    uint_t cache_depth;
    gin_bs_read_word(bs, ridx, GIN_GIN_CACHE_DEPTH_BIT_LENGTH, &cache_depth);
    cache->depth = (int_t)cache_depth;
    ridx += GIN_GIN_CACHE_DEPTH_BIT_LENGTH;
    // read the number of entries
    uint_t no_entries;
    gin_bs_read_word(bs, ridx, GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH, &no_entries);
    cache->no_entries = (int_t)no_entries;
    ridx += GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH;
    // fmi size in bytes
    uint_t key_fmi_size_in_bits;
    gin_bs_read_word(bs, ridx, GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH, &key_fmi_size_in_bits);
    cache->key_fmi_size_in_bits = (int_t)key_fmi_size_in_bits;
    ridx += GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH;
    // value buffer size in bytes
    uint_t value_buffer_size_in_bits;
    gin_bs_read_word(bs, ridx, GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH, &value_buffer_size_in_bits);
    cache->value_buffer_size_in_bits = (int_t)value_buffer_size_in_bits;
    ridx += GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH;
    // word align
    ridx  = (1 + ((ridx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 1 - Parse item word offsets
    **************************************************************************/
    cache->item_offsets = calloc(cache->no_entries, sizeof(word_t));
    memcpy(cache->item_offsets, bs->words + (ridx >> WORD_LOG_BITS), cache->no_entries * sizeof(word_t));
    ridx += no_entries * WORD_NUM_BITS;
    // word align
    ridx  = (1 + ((ridx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 2 - Parse FMI bitstream
    **************************************************************************/

    word_t fmi_no_bits = 0;
    gin_bs_read_word(bs, ridx, GIN_GIN_FMI_NO_BITS_BIT_LENGTH, &fmi_no_bits);
    ridx += GIN_GIN_FMI_NO_BITS_BIT_LENGTH;
    // read padded words into buffer
    int_t fmi_buf_size_bytes = (int_t)fmi_no_bits >> 3;
    int_t fmi_buf_size_words_padded = 1 + (((int_t)fmi_no_bits - 1) >> WORD_LOG_BITS);
    uint8_t *fmi_buf = calloc(fmi_buf_size_words_padded, sizeof(word_t));
    uint64_t fmi_buf_size = fmi_buf_size_bytes;
    word_t *fmi_buf_word = (word_t*)fmi_buf;
    for(int_t i = 0; i < fmi_buf_size_words_padded; i++) {
        word_t word = 0;
        gin_bs_read_word(bs, ridx, WORD_NUM_BITS, &word);
        fmi_buf_word[i] = word;
        ridx += WORD_NUM_BITS;
    }
#ifdef GIN_SDSL
    cache->key_fmi = csa_wt_from_buffer(fmi_buf, fmi_buf_size);
#else
    gin_fmi_serialize_from_buffer(fmi_buf, fmi_buf_size, &cache->key_fmi);
#endif
    free(fmi_buf);
    // word align
    ridx  = (1 + ((ridx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;

    /******************************************************
    * Step 3 - Parse the values
    ******************************************************/
    word_t no_value_buffer_words = (word_t)(value_buffer_size_in_bits) >> WORD_LOG_BITS;
    cache->items = calloc(no_value_buffer_words, sizeof(word_t));
    memcpy(cache->items, bs->words + (ridx >> WORD_LOG_BITS), no_value_buffer_words * sizeof(word_t));
    /******************************************************
    * Step 4 - Cleanup
    ******************************************************/
    gin_bs_free(bs);
    *cachew = cache;
}

void gin_gin_cache_free(gin_gin_cache_t *cache) {
    if(!cache) return;
#ifdef GIN_SDSL
    csa_wt_free(cache->key_fmi);
#else
    gin_fmi_free(cache->key_fmi);
#endif
    free(cache->items);
    free(cache->item_offsets);
    free(cache);
}

bool gin_gin_advance_fork(gin_gin_t *gin, gin_fork_node_t *fork, gin_string_t *pattern) {
#ifdef GIN_SDSL
    sdsl_csa *fmi = gin->graph_fmi;
    char c = pattern->seq[fork->pos];
    count_t rank_lo = fork->sa_lo ? csa_wt_rank(fmi, fork->sa_lo, c) : 0ull;
    count_t rank_hi = csa_wt_rank(fmi, fork->sa_hi, c);
    uint64_t base = csa_wt_char_sa_base(fmi, c);
    fork->sa_lo = (int_t)(base + rank_lo);
    fork->sa_hi = (int_t)(base + rank_hi);
    --fork->pos;
    return fork->sa_hi > fork->sa_lo;
#else
    gin_fmi_t *fmi = gin->graph_fmi;
    word_t encoding = fmi->c2e[pattern->seq[fork->pos]];
    count_t rank_lo_m_1 = fork->sa_lo ? gin_fmi_rank(fmi,encoding, fork->sa_lo-1) : 0ull;
    count_t rank_hi_m_1 = fork->sa_hi ? gin_fmi_rank(fmi,encoding, fork->sa_hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    fork->sa_lo = (int_t)(base + rank_lo_m_1);
    fork->sa_hi = (int_t)(base + rank_hi_m_1);
    --fork->pos;
    return fork->sa_hi > fork->sa_lo;
#endif
}

bool gin_gin_fork_precedence_range(gin_gin_t *gin, gin_fork_node_t *fork, char_t c, int_t *lo, int_t *hi) {
#ifdef GIN_SDSL
    sdsl_csa *fmi = gin->graph_fmi;
    count_t rank_lo = fork->sa_lo ? csa_wt_rank(fmi, fork->sa_lo, c) : 0ull;
    count_t rank_hi = csa_wt_rank(fmi, fork->sa_hi, c);
    uint64_t base = csa_wt_char_sa_base(fmi, c);
    *lo = (int_t)(base + rank_lo);
    *hi = (int_t)(base + rank_hi);
    return *hi > *lo;
#else
    gin_fmi_t *fmi = gin->graph_fmi;
    word_t encoding = fmi->c2e[c];
    count_t rank_lo_m_1 = fork->sa_lo ? gin_fmi_rank(fmi,encoding, fork->sa_lo-1) : 0ull;
    count_t rank_hi_m_1 = fork->sa_hi ? gin_fmi_rank(fmi,encoding, fork->sa_hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    *lo = (int_t)(base + rank_lo_m_1);
    *hi = (int_t)(base + rank_hi_m_1);
    return *hi > *lo;
#endif
}

gin_vector_t *gin_gin_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords) {
    gin_vector_t *rval;
    gin_vector_init(&rval,no_codewords,&gin_fstruct_string);
    int_t len = (int_t)gin_ceil_log2(no_codewords); // length of each codeword
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
        gin_string_t *copy;
        gin_string_init_cstr(&copy, codeword);
        gin_vector_append(rval, copy);
    }
    return rval;
}

void gin_gin_serialize_from_buffer(gin_gin_t **gin_ret, unsigned char *buf, uint64_t buf_size) {
    gin_gin_t *gin = calloc(1, sizeof(gin_gin_t));
    if(!gin) {
        *gin_ret = NULL;
        return;
    }
    /**************************************************************************
    * Step 0 - Wrap the buffer into a bitstream - this consumes the buffer
    **************************************************************************/
    gin_bs_t *bs;
    gin_bs_init_from_buffer(buf, buf_size, &bs);
    uint_t ridx = 0;
    /**************************************************************************
    * Step 1 - Read total size and special characters
    **************************************************************************/
    word_t gin_size_in_bits, c_0, c_1;
    gin_bs_read_word(bs, ridx, GIN_GIN_NO_BITS_BIT_LENGTH, &gin_size_in_bits);
    ridx += GIN_GIN_NO_BITS_BIT_LENGTH;
    gin_bs_read_word(bs, ridx, GIN_GIN_SPECIAL_CHAR_BIT_LENGTH, &c_0);
    ridx += GIN_GIN_SPECIAL_CHAR_BIT_LENGTH;
    gin_bs_read_word(bs, ridx, GIN_GIN_SPECIAL_CHAR_BIT_LENGTH, &c_1);
    ridx += GIN_GIN_SPECIAL_CHAR_BIT_LENGTH;
    gin->c_0 = (char_t)c_0;
    gin->c_1 = (char_t)c_1;
    /**************************************************************************
    * Step 2 - Read the permutation and BWT to vertex ID map
    **************************************************************************/
    word_t no_vertices;
    gin_bs_read_word(bs, ridx, GIN_GIN_NO_VERTICES_BIT_LENGTH, &no_vertices);
    ridx += GIN_GIN_NO_VERTICES_BIT_LENGTH;

    gin_vector_t *permutation, *bwt_to_vid;
    gin_vector_init(&permutation, (int_t)no_vertices, &prm_fstruct);
    gin_vector_init(&bwt_to_vid, (int_t)no_vertices, &prm_fstruct);
    word_t p;
    for(int_t i = 0; i < no_vertices; i++) {
        gin_bs_read_word(bs, ridx, GIN_GIN_PERMUTATION_BIT_LENGTH, &p);
        ridx += GIN_GIN_PERMUTATION_BIT_LENGTH;
        gin_vector_append(permutation, (void*)p);
    }
    for(int_t i = 0; i < no_vertices; i++) {
        gin_bs_read_word(bs, ridx, GIN_GIN_BWT_TO_VID_BIT_LENGTH, &p);
        ridx += GIN_GIN_BWT_TO_VID_BIT_LENGTH;
        gin_vector_append(bwt_to_vid, (void*)p);
    }
    gin_vector_fit(permutation);
    gin_vector_fit(bwt_to_vid);
    gin->permutation = permutation;
    gin->bwt_to_vid = bwt_to_vid;
    /**************************************************************************
    * Step 3 - Reconstruct the FMI from the bitstream
    **************************************************************************/

    word_t fmi_no_bits = 0;
    gin_bs_read_word(bs, ridx, GIN_GIN_FMI_NO_BITS_BIT_LENGTH, &fmi_no_bits);
    ridx += GIN_GIN_FMI_NO_BITS_BIT_LENGTH;
    // read padded words into buffer
    int_t fmi_buf_size_bytes = (int_t)fmi_no_bits >> 3;
    int_t fmi_buf_size_words_padded = 1 + (((int_t)fmi_no_bits - 1) >> WORD_LOG_BITS);
    uint8_t *fmi_buf = calloc(fmi_buf_size_words_padded, sizeof(word_t));
    uint64_t fmi_buf_size = fmi_buf_size_bytes;
    word_t *fmi_buf_word = (word_t*)fmi_buf;
    for(int_t i = 0; i < fmi_buf_size_words_padded; i++) {
        word_t word = 0;
        gin_bs_read_word(bs, ridx, WORD_NUM_BITS, &word);
        fmi_buf_word[i] = word;
        ridx += WORD_NUM_BITS;
    }
#ifdef GIN_SDSL
    gin->graph_fmi = csa_wt_from_buffer(fmi_buf, fmi_buf_size);
    gin->no_chars = csa_wt_bwt_length(gin->graph_fmi);
    csa_wt_populate_alphabet(gin->graph_fmi, &gin->alphabet, &gin->alphabet_size);
#else
    gin_fmi_serialize_from_buffer(fmi_buf, fmi_buf_size, &gin->graph_fmi);
    gin->no_chars = gin->graph_fmi->no_chars;
    gin->alphabet_size = gin->graph_fmi->alphabet_size;
    gin->alphabet = calloc(gin->alphabet_size, sizeof(int_t));
    memcpy(gin->alphabet, gin->graph_fmi->alphabet, sizeof(int_t) * gin->alphabet_size);
#endif
    free(fmi_buf);
    // word align
    ridx  = (1 + ((ridx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 4 - Construct the interval-merge tree from the leaves
    **************************************************************************/
    word_t gin_imt_bit_length;
    gin_bs_read_word(bs, ridx, GIN_GIN_IMT_NO_BITS_BIT_LENGTH, &gin_imt_bit_length);
    ridx += GIN_GIN_IMT_NO_BITS_BIT_LENGTH;
    // no keys is equal to the number of vertices
    gin_vector_t *kv_pairs;
    gin_vector_init(&kv_pairs, (int_t)no_vertices, &gin_fstruct_vector);
    word_t no_intervals, lo, hi;
    for(int_t i = 0; i < no_vertices; i++) {
        gin_bs_read_word(bs, ridx, GIN_GIN_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH, &no_intervals);
        ridx += GIN_GIN_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH;
        gin_vector_t *bwt_neighbor_intervals;
        gin_vector_init(&bwt_neighbor_intervals, (int_t)no_intervals, &gin_fstruct_imt_interval);
        for(int_t j = 0; j < (int_t)no_intervals; j++) {
            gin_bs_read_word(bs, ridx, GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH, &lo);
            ridx += GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            gin_bs_read_word(bs, ridx, GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH, &hi);
            ridx += GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            gin_imt_interval_t *interval;
            gin_imt_interval_init(&interval, (int_t)lo, (int_t)hi);
            gin_vector_append(bwt_neighbor_intervals, interval);
        }
        gin_vector_append(kv_pairs, bwt_neighbor_intervals);
    }
    gin_imt_t *inverval_merge_tree;
    gin_imt_init(&inverval_merge_tree, (int_t)no_vertices, kv_pairs);
    gin->r2r_tree = inverval_merge_tree;
    #ifdef GIN_ORACLE
    int_t V = gin->permutation->size;
    int_t *alphabet = gin->alphabet;
    int_t alphabet_size = gin->alphabet_size;
    int_t *vertex_last_char_enc = calloc(V, sizeof(int_t));
    #ifdef GIN_SDSL
    csa_wt_bwt(gin->graph_fmi, (uint64_t*)vertex_last_char_enc, V+1, 2*V);
    #else
    gin_fmi_bwt(gin->graph_fmi, (uint64_t*)vertex_last_char_enc, V+1, 2*V);
    #endif
    int_t *c2e = calloc(GIN_MAX_ALPHABET_SIZE, sizeof(int_t));
    int_t t = 0;
    for(int_t i = 0; i < alphabet_size; i++) {
        c2e[alphabet[i]] = t;
        t++;
    }
    for(int_t i = 0; i < V; i++) {
        vertex_last_char_enc[i] = c2e[vertex_last_char_enc[i]];
    }
    free(c2e);
    gin_oimt_init(gin->r2r_tree, vertex_last_char_enc, alphabet, (int_t)alphabet_size, &gin->oracle_r2r);
    free(vertex_last_char_enc);
    #endif
    gin_vector_free(kv_pairs);
    /**************************************************************************
    * Step 5 - Cleanup and return the reconstructed index
    **************************************************************************/
    gin_bs_free(bs);
    *gin_ret = gin;
}

void gin_gin_serialize_to_buffer(gin_gin_t *gin, unsigned char **buf_ret, uint64_t *buf_size_ret) {
    if(!gin) {
        *buf_ret = NULL;
        *buf_size_ret = 0;
        return;
    }
    // init bitstream to write everything
    gin_bs_t *bs;
    gin_bs_init(&bs);
    uint_t widx = 0;
    /**************************************************************************
    * Step 0 - Write the special characters to the buffer
    **************************************************************************/
    /* Reserve the beginning of the buffer for the total bit length */
    widx += GIN_GIN_NO_BITS_BIT_LENGTH;
    gin_bs_write_word(bs, widx, gin->c_0, GIN_GIN_SPECIAL_CHAR_BIT_LENGTH);
    widx += GIN_GIN_SPECIAL_CHAR_BIT_LENGTH;
    gin_bs_write_word(bs, widx, gin->c_1, GIN_GIN_SPECIAL_CHAR_BIT_LENGTH);
    widx += GIN_GIN_SPECIAL_CHAR_BIT_LENGTH;
    /**************************************************************************
    * Step 1 - Write the permutation to the buffer
    **************************************************************************/
    gin_bs_write_word(bs, widx, gin->permutation->size, GIN_GIN_NO_VERTICES_BIT_LENGTH);
    widx += GIN_GIN_NO_VERTICES_BIT_LENGTH;
    for(int_t i = 0; i < gin->permutation->size; i++) {
        gin_bs_write_word(bs, widx, (word_t)gin->permutation->data[i], GIN_GIN_PERMUTATION_BIT_LENGTH);
        widx += GIN_GIN_PERMUTATION_BIT_LENGTH;
    }
    /**************************************************************************
    * Step 2 - Write the BWT c_0 to VID mapping
    **************************************************************************/
    for(int_t i = 0; i < gin->bwt_to_vid->size; i++) {
        gin_bs_write_word(bs, widx, (word_t)gin->bwt_to_vid->data[i], GIN_GIN_BWT_TO_VID_BIT_LENGTH);
        widx += GIN_GIN_BWT_TO_VID_BIT_LENGTH;
    }
    /**************************************************************************
    * Step 3 - Write the FMI of the graph encoding
    **************************************************************************/

    uint8_t *fmi_buf = NULL;
    uint64_t fmi_buf_size = 0;
#ifdef GIN_SDSL
    csa_wt_to_buffer(gin->graph_fmi, &fmi_buf, &fmi_buf_size);
#else
    gin_fmi_serialize_to_buffer(gin->graph_fmi, &fmi_buf, &fmi_buf_size);
#endif
    int_t fmi_buf_size_in_bits = (int_t)fmi_buf_size << 3;
    gin_bs_write_word(bs, widx, (word_t)fmi_buf_size_in_bits, GIN_GIN_FMI_NO_BITS_BIT_LENGTH);
    widx += GIN_GIN_FMI_NO_BITS_BIT_LENGTH;
    int_t padded_no_bits = (1 + ((fmi_buf_size_in_bits-1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    int_t fmi_buf_no_padded_words = padded_no_bits >> WORD_LOG_BITS;
    word_t *fmi_buf_word = (word_t*)fmi_buf;
    for(int_t i = 0; i < fmi_buf_no_padded_words; i++) {
        gin_bs_write_word(bs, widx, fmi_buf_word[i], WORD_NUM_BITS);
        widx += WORD_NUM_BITS;
    }
    free(fmi_buf);
    // word align
    widx  = (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /**************************************************************************
    * Step 4 - Write the IMT bitstream
    **************************************************************************/
    /* Reserve space for the IMT encoding bit length */
    widx += GIN_GIN_IMT_NO_BITS_BIT_LENGTH;
    uint_t prev_widx = widx;
    gin_gin_serialize_to_buffer_imt_helper(gin->r2r_tree->root, bs, &widx);
    gin_bs_write_word(bs, prev_widx - GIN_GIN_IMT_NO_BITS_BIT_LENGTH, widx - prev_widx, GIN_GIN_IMT_NO_BITS_BIT_LENGTH);
    /**************************************************************************
    * Step 5 - Write total length in the beginning
    **************************************************************************/
    gin_bs_write_word(bs, 0, widx, GIN_GIN_NO_BITS_BIT_LENGTH);
    /**************************************************************************
    * Step 6 - Return
    **************************************************************************/
    gin_bs_fit(bs, widx);
    word_t *buf;
    uint_t no_words;
    gin_bs_detach(bs, &buf, &no_words);
    gin_bs_free(bs);
    *buf_ret = (unsigned char*)buf;
    *buf_size_ret = no_words * sizeof(word_t);
}

void gin_gin_serialize_to_buffer_imt_helper(gin_imt_node_t *node, gin_bs_t *bs, uint_t *widx) {
    if(node) {
        gin_gin_serialize_to_buffer_imt_helper((gin_imt_node_t*)node->left, bs, widx);
        if (node->hi == node->lo) {
            gin_bs_write_word(bs, *widx, (word_t)node->intervals->size, GIN_GIN_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH);
            *widx += GIN_GIN_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH;
            for(int_t i = 0; i < node->intervals->size; i++) {
                int_t lo = ((gin_imt_interval_t*)(node->intervals->data[i]))->lo;
                int_t hi = ((gin_imt_interval_t*)(node->intervals->data[i]))->hi;
                gin_bs_write_word(bs, *widx, lo, GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH);
                *widx += GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
                gin_bs_write_word(bs, *widx, hi, GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH);
                *widx += GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH;
            }
        }
        gin_gin_serialize_to_buffer_imt_helper((gin_imt_node_t*)node->right, bs, widx);
    }
}