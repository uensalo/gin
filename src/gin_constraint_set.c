/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_constraint_set.c is part of gin
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
#include "gin_constraint_set.h"

void gin_constraint_set_enumerate(gin_vector_t **constraint_sets, gin_graph_t *graph, int_t depth, bool multiple_vertex_span) {
    gin_table_t *c_table = gin_constraint_set_extract(graph, depth, multiple_vertex_span);
    gin_vector_t *constraints = gin_constraint_set_to_vector(c_table);
    // soft-free the table
    c_table->key_f = &prm_fstruct;
    c_table->val_f = &prm_fstruct;
    gin_table_free(c_table);
    // return the sorted constraints
    *constraint_sets = constraints;
}

void gin_constraint_set_free(gin_constraint_set_t *cs) {
    free(cs);
}
uint_t gin_constraint_set_hash(gin_constraint_set_t *c1) {
    return gin_string_hash(c1->str);
}
int gin_constraint_set_comp(gin_constraint_set_t *c1, gin_constraint_set_t *c2) {
    if (c1->str->size < c2->str->size) {
        return -1;
    } else if (c1->str->size > c2->str->size){
        return 1;
    }
    return gin_string_comp(c1->str, c2->str);
}
void *gin_constraint_set_copy(gin_constraint_set_t *c1) {
    gin_constraint_set_t *copy = calloc(1, sizeof(gin_constraint_set_t));
    copy->str = gin_string_copy(c1->str);
    copy->vertices = gin_vector_copy(c1->vertices);
    return copy;
}

void gin_constraint_set_flatten_alphabet(void *key, void *value, void *params) {
    gin_graph_flatten_alphabet_t *p = (gin_graph_flatten_alphabet_t*)params;
    int_t idx = p->alphabet->size;
    gin_table_insert(p->char2idx, key, (void*)idx);
    gin_table_insert(p->idx2char, (void*)idx, key);
    gin_vector_append(p->alphabet, key);
}

gin_table_t *gin_constraint_set_extract(gin_graph_t *graph, int_t max_depth, bool multiple_vertex_span) {
    // first, construct an alphabet and encoding translation tables
    gin_vector_t *alphabet;
    gin_tree_t *alphabet_tree;
    gin_table_t *char2idx, *idx2char;
    gin_tree_init(&alphabet_tree, &prm_fstruct, &prm_fstruct);
    // traverse all vertices and add characters to a set
    for(int_t i = 0; i < graph->vertex_list->size; i++) {
        void *vertex;
        gin_table_lookup(graph->vertices, (void*)i, &vertex);
        gin_string_t *label = ((gin_vertex_t*)vertex)->label;
        for(int j = 0; j < label->size; j++) {
            gin_tree_insert(alphabet_tree, (void*)label->seq[j], NULL);
        }
    }
    // flatten the tree set to obtain the alphabet and translation tables
    gin_graph_flatten_alphabet_t p;
    gin_vector_init(&alphabet, GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    gin_table_init(&char2idx, GIN_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    gin_table_init(&idx2char, GIN_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    p.alphabet = alphabet;
    p.char2idx = char2idx;
    p.idx2char = idx2char;
    gin_tree_inorder(alphabet_tree, &p, gin_constraint_set_flatten_alphabet);
    gin_tree_free(alphabet_tree);

    // initialize the recursion
    gin_string_t *initial_prefix;
    gin_string_init(&initial_prefix, GIN_STRING_INIT_SIZE);

    gin_table_t *constraint_sets;
    gin_table_init(&constraint_sets, GIN_HT_INIT_SIZE, &gin_fstruct_string, &gin_fstruct_vector);

    gin_vector_t *paths;
    gin_vector_init(&paths, graph->vertices->size, &prm_fstruct);
    for(int_t i = 0; i < graph->vertices->size; i++) {
        gin_graph_ecs_t *path = calloc(1, sizeof(gin_graph_ecs_t));
        path->pos = 0;
        path->head_vid = (vid_t)i;
        path->end_vid = (vid_t)i;
        gin_vector_append(paths, path);
    }
    gin_constraint_set_extract_helper(paths, initial_prefix, constraint_sets, alphabet, char2idx, idx2char, graph, max_depth, multiple_vertex_span);

    // cleanup and return
    gin_table_free(char2idx);
    gin_table_free(idx2char);
    gin_vector_free(alphabet);
    gin_string_free(initial_prefix); // free the initial prefix
    return constraint_sets;
}

void gin_constraint_set_extract_helper(gin_vector_t *paths,
                                       gin_string_t *prefix,
                                       gin_table_t *constraint_sets,
                                       gin_vector_t *alphabet,
                                       gin_table_t *char2idx,
                                       gin_table_t *idx2char,
                                       gin_graph_t *graph,
                                       int_t max_depth,
                                       bool multiple_vertex_span) {
    if(!paths->size) {
        // empty bucket - free the relevant stuff and return without doing anything
        gin_vector_free(paths);
        gin_string_free(prefix);
        return;
    }
    // sort paths into buckets
    gin_vector_t **buckets = calloc(alphabet->size, sizeof(gin_vector_t*));
    // init buckets
    for(int_t i = 0; i < alphabet->size; i++) {
        gin_vector_init(&buckets[i], GIN_VECTOR_INIT_SIZE, &prm_fstruct);
    }

    for(int_t i = 0; i < paths->size; i++) {
        gin_graph_ecs_t *path = paths->data[i];
        void *last_vertex;
        gin_table_lookup(graph->vertices, (void*)path->end_vid, &last_vertex);
        gin_string_t *label = ((gin_vertex_t*)last_vertex)->label;
        if(multiple_vertex_span && path->pos >= label->size) {
            // path exhausted! Consider all outgoing neighbors and generate new paths from this path...
            gin_vector_t *outgoing_neighbors;
            gin_table_lookup(graph->outgoing_neighbors, (void*)path->end_vid, &outgoing_neighbors);
            for(int_t j = 0; j < outgoing_neighbors->size; j++) {
                gin_graph_ecs_t *new_path = calloc(1, sizeof(gin_graph_ecs_t));
                new_path->head_vid = path->head_vid;
                new_path->end_vid = (vid_t)outgoing_neighbors->data[j];
                new_path->pos = 0;
                // put the path into its respective bucket
                void *new_vertex;
                gin_table_lookup(graph->vertices, (void*)new_path->end_vid, &new_vertex);
                char_t last_char = ((gin_vertex_t*)new_vertex)->label->seq[new_path->pos];
                void *bucket_idx;
                gin_table_lookup(char2idx, (void*)last_char, &bucket_idx);
                gin_vector_append(buckets[(int_t)bucket_idx], new_path);
                ++new_path->pos;
            }
            // the old path is no longer useful as it will not be added anywhere, kill it
            free(path);
        } else {
            char_t last_char = label->seq[path->pos];
            void *bucket_idx;
            gin_table_lookup(char2idx, (void*)last_char, &bucket_idx);
            gin_vector_append(buckets[(int_t)bucket_idx], path);
            ++path->pos;
        }
    }
    // extract constraint sets and call recursively
    for(int_t i = 0; i < alphabet->size; i++) {
        if(!buckets[i]->size) continue;
        // compute bucket prefix
        gin_string_t *bucket_prefix = gin_string_copy(prefix);
        void* bucket_char;
        gin_table_lookup(idx2char, (void*)i, &bucket_char);
        gin_string_append(bucket_prefix, (char_t)bucket_char);

        // find the union of all path roots of the bucket
        gin_tree_t *constraint_values;
        gin_tree_init(&constraint_values, &prm_fstruct, &prm_fstruct);
        for(int j = 0; j < buckets[i]->size; j++) {
            gin_graph_ecs_t *path = buckets[i]->data[j];
            gin_vector_t *incoming_neighbors;
            gin_table_lookup(graph->incoming_neighbors, (void*)path->head_vid, &incoming_neighbors);
            for(int k = 0; k < incoming_neighbors->size; k++) {
                gin_tree_insert(constraint_values, incoming_neighbors->data[k], incoming_neighbors->data[k]);
            }
        }
        // flatten keys of the tree
        gin_vector_t *constraints;
        gin_vector_init(&constraints, constraint_values->no_items, &prm_fstruct);
        gin_tree_inorder(constraint_values, constraints, gin_constraint_set_flatten);

        // add constraint set to constraint set table, with bucket prefix as key
        gin_table_insert(constraint_sets, bucket_prefix, constraints);
        gin_tree_free(constraint_values);

        // call to children recursively
        if(bucket_prefix->size < max_depth) {
            gin_constraint_set_extract_helper(buckets[i], bucket_prefix, constraint_sets, alphabet, char2idx,
                                              idx2char, graph, max_depth, multiple_vertex_span);
        }
    }
    // cleanup buckets and paths
    for(int_t i = 0; i < alphabet->size; i++) {
        for(int_t j = 0; j < buckets[i]->size; j++) {
    //        free(buckets[i]->data[j]);
        }
        gin_vector_free(buckets[i]);
    }
    free(buckets);
}

void gin_constraint_set_flatten(void *key, void *value, void *p) {
    gin_vector_t *constraints = (gin_vector_t*)p;
    gin_vector_append(constraints, key);
}

void gin_constraint_set_to_vector_helper(void *key, void* value, void *p) {
    gin_vector_t *constraints = (gin_vector_t*)p;
    gin_constraint_set_t *constraint = calloc(1, sizeof(gin_constraint_set_t));
    constraint->str = key;
    constraint->vertices = value;
    gin_vector_append(constraints, constraint);
}

gin_vector_t *gin_constraint_set_to_vector(gin_table_t *constraint_sets) {
    gin_vector_t *constraints;
    gin_vector_init(&constraints, GIN_VECTOR_INIT_SIZE, &gin_fstruct_cs);
    gin_table_traverse(constraint_sets, constraints, gin_constraint_set_to_vector_helper);
    gin_vector_sort(constraints);
    return constraints;
}