#include "fmd_graph.h"

void fmd_vertex_init(fmd_vertex_t **vertex, vid_t id, fmd_string_t *label) {
    fmd_vertex_t *v = calloc(1, sizeof(fmd_vertex_t));
    if(!v) {
        *vertex = NULL;
        return;
    }
    v->id = id;
    v->label = label;
    *vertex = v;
}

void fmd_vertex_free(fmd_vertex_t *vertex) {
    fmd_string_free(vertex->label);
    free(vertex);
}

uint_t fmd_vertex_hash(fmd_vertex_t *vertex) {
    // djb2
    uint_t hash = 5381;
    hash = ((hash << 5) + hash) + (size_t)vertex->id;
    uint_t label_hash = fmd_string_hash(vertex->label);
    hash = ((hash << 5) + hash) + label_hash;
    return hash;
}

int fmd_vertex_comp(fmd_vertex_t *v1, fmd_vertex_t *v2) {
    if (v1->id != v2->id) {
        return v1->id < v2->id ? -1 : 1;
    }
    return fmd_string_comp(v1->label, v2->label);
}

fmd_vertex_t *fmd_vertex_copy(fmd_vertex_t *vertex) {
    fmd_vertex_t *copy = calloc(1,sizeof(fmd_vertex_t));
    if(!copy) return NULL;
    copy->id = vertex->id;
    copy->label = fmd_string_copy(vertex->label);
    return copy;
}

void fmd_graph_init(fmd_graph_t **graph) {
    *graph = calloc(1, sizeof(fmd_graph_t));
    fmd_vector_init(&(*graph)->vertex_list, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_init(&(*graph)->vertices, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vertex);
    fmd_table_init(&(*graph)->incoming_neighbors, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vector);
    fmd_table_init(&(*graph)->outgoing_neighbors, FMD_HT_INIT_SIZE, &prm_fstruct, &fmd_fstruct_vector);
}

void fmd_graph_free(fmd_graph_t *graph) {
    fmd_vector_free(graph->vertex_list);
    fmd_table_free(graph->vertices);
    fmd_table_free(graph->incoming_neighbors);
    fmd_table_free(graph->outgoing_neighbors);
    free(graph);
}

void fmd_graph_insert_vertex(fmd_graph_t *graph, vid_t id, fmd_string_t *label) {
    fmd_vertex_t *vertex;
    fmd_vertex_init(&vertex, id, label);
    fmd_table_insert(graph->vertices, (void*)id, vertex);
    fmd_vector_append(graph->vertex_list, vertex);

    fmd_vector_t *incoming_neighbors;
    fmd_vector_init(&incoming_neighbors, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_insert(graph->incoming_neighbors, (void*)id, incoming_neighbors);

    fmd_vector_t *outgoing_neighbors;
    fmd_vector_init(&outgoing_neighbors, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_insert(graph->outgoing_neighbors, (void*)id, outgoing_neighbors);
}

void fmd_graph_insert_edge(fmd_graph_t *graph, vid_t source, vid_t destination) {
    void *incoming_neighbors;
    fmd_table_lookup(graph->incoming_neighbors, (void*)destination, &incoming_neighbors);
    fmd_vector_append((fmd_vector_t*)incoming_neighbors, (void*)source);

    void *outgoing_neigbors;
    fmd_table_lookup(graph->outgoing_neighbors, (void*)source, &outgoing_neigbors);
    fmd_vector_append((fmd_vector_t*)outgoing_neigbors, (void*)destination);

    ++graph->no_edges;
}

uint_t fmd_graph_hash(fmd_graph_t *graph) {
    uint_t vertices_hash = fmd_table_hash(graph->vertices);
    uint_t edges_hash = fmd_table_hash(graph->incoming_neighbors);
    return (vertices_hash * 31) + edges_hash;
}

int fmd_graph_comp(fmd_graph_t *g1, fmd_graph_t *g2) {
    int vertices_comp = fmd_table_comp(g1->vertices, g2->vertices);
    if (vertices_comp != 0) {
        return vertices_comp;
    }
    return fmd_table_comp(g1->incoming_neighbors, g2->incoming_neighbors);
}

fmd_graph_t *fmd_graph_copy(fmd_graph_t *graph) {
    fmd_graph_t *copy = calloc(1, sizeof(fmd_graph_t));
    copy->vertex_list = fmd_vector_copy(graph->vertex_list);
    copy->vertices = fmd_table_copy(graph->vertices);
    copy->incoming_neighbors = fmd_table_copy(graph->incoming_neighbors);
    copy->outgoing_neighbors = fmd_table_copy(graph->outgoing_neighbors);
    return copy;
}

void fmd_graph_flatten_alphabet(void *key, void *value, void *params) {
    fmd_graph_flatten_alphabet_t *p = (fmd_graph_flatten_alphabet_t*)params;
    int_t idx = p->alphabet->size;
    fmd_table_insert(p->char2idx, key, (void*)idx);
    fmd_table_insert(p->idx2char, (void*)idx, key);
    fmd_vector_append(p->alphabet, key);
}

fmd_table_t *fmd_graph_extract_constraint_sets(fmd_graph_t *graph, int_t max_depth) {
    // first, construct an alphabet and encoding translation tables
    fmd_vector_t *alphabet;
    fmd_tree_t *alphabet_tree;
    fmd_table_t *char2idx, *idx2char;
    fmd_tree_init(&alphabet_tree, &prm_fstruct, &prm_fstruct);
    // traverse all vertices and add characters to a set
    for(int_t i = 0; i < graph->vertex_list->size; i++) {
        void *vertex;
        fmd_table_lookup(graph->vertices, (void*)i, &vertex);
        fmd_string_t *label = ((fmd_vertex_t*)vertex)->label;
        for(int j = 0; j < label->size; j++) {
            fmd_tree_insert(alphabet_tree, (void*)label->seq[j], NULL);
        }
    }
    // flatten the tree set to obtain the alphabet and translation tables
    fmd_graph_flatten_alphabet_t p;
    fmd_vector_init(&alphabet, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    fmd_table_init(&char2idx, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    fmd_table_init(&idx2char, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    p.alphabet = alphabet;
    p.char2idx = char2idx;
    p.idx2char = idx2char;
    fmd_tree_inorder(alphabet_tree, &p, fmd_graph_flatten_alphabet);
    fmd_tree_free(alphabet_tree);

    // initialize the recursion
    fmd_string_t *initial_prefix;
    fmd_string_init(&initial_prefix, FMD_STRING_INIT_SIZE);

    fmd_table_t *constraint_sets;
    fmd_table_init(&constraint_sets, FMD_HT_INIT_SIZE, &fmd_fstruct_string, &fmd_fstruct_vector);

    fmd_vector_t *paths;
    fmd_vector_init(&paths, graph->vertices->size, &prm_fstruct);
    for(int_t i = 0; i < graph->vertices->size; i++) {
        fmd_graph_ecs_t *path = calloc(1, sizeof(fmd_graph_ecs_t));
        path->pos = 0;
        path->head_vid = (vid_t)i;
        path->end_vid = (vid_t)i;
        fmd_vector_append(paths, path);
    }
    fmd_graph_extract_constraint_sets_helper(paths, initial_prefix, constraint_sets, alphabet, char2idx, idx2char, graph, max_depth);

    // cleanup and return
    fmd_table_free(char2idx);
    fmd_table_free(idx2char);
    fmd_vector_free(alphabet);
    fmd_string_free(initial_prefix); // free the initial prefix
    return constraint_sets;
}



void fmd_graph_extract_constraint_sets_helper(fmd_vector_t *paths,
                                              fmd_string_t *prefix,
                                              fmd_table_t *constraint_sets,
                                              fmd_vector_t *alphabet,
                                              fmd_table_t *char2idx,
                                              fmd_table_t *idx2char,
                                              fmd_graph_t *graph,
                                              int_t max_depth) {
    // sort paths into buckets
    fmd_vector_t **buckets = calloc(alphabet->size, sizeof(fmd_vector_t*));
    // init buckets
    for(int_t i = 0; i < alphabet->size; i++) {
        fmd_vector_init(&buckets[i], FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    }

    for(int_t i = 0; i < paths->size; i++) {
        fmd_graph_ecs_t *path = paths->data[i];
        void *last_vertex;
        fmd_table_lookup(graph->vertices, (void*)path->end_vid, &last_vertex);
        fmd_string_t *label = ((fmd_vertex_t*)last_vertex)->label;
        if(path->pos >= label->size) {
            // path exhausted! Consider all outgoing neighbors and generate new paths from this path...
            fmd_vector_t *outgoing_neighbors;
            fmd_table_lookup(graph->outgoing_neighbors, (void*)path->end_vid, &outgoing_neighbors);
            for(int_t j = 0; j < outgoing_neighbors->size; j++) {
                fmd_graph_ecs_t *new_path = calloc(1, sizeof(fmd_graph_ecs_t));
                new_path->head_vid = path->head_vid;
                new_path->end_vid = (vid_t)outgoing_neighbors->data[j];
                new_path->pos = 0;
                // put the path into its respective bucket
                void *new_vertex;
                fmd_table_lookup(graph->vertices, (void*)new_path->end_vid, &new_vertex);
                char_t last_char = ((fmd_vertex_t*)new_vertex)->label->seq[new_path->pos];
                void *bucket_idx;
                fmd_table_lookup(char2idx, (void*)last_char, &bucket_idx);
                fmd_vector_append(buckets[(int_t)bucket_idx], new_path);
                ++new_path->pos;
            }
            // the old path is no longer useful as it will not be added anywhere, kill it
            free(path);
        } else {
            char_t last_char = label->seq[path->pos];
            void *bucket_idx;
            fmd_table_lookup(char2idx, (void*)last_char, &bucket_idx);
            fmd_vector_append(buckets[(int_t)bucket_idx], path);
            ++path->pos;
        }
    }
    // extract constraint sets and call recursively
    for(int_t i = 0; i < alphabet->size; i++) {
        if(!buckets[i]->size) continue;
        // compute bucket prefix
        fmd_string_t *bucket_prefix = fmd_string_copy(prefix);
        void* bucket_char;
        fmd_table_lookup(idx2char, (void*)i, &bucket_char);
        fmd_string_append(bucket_prefix, (char_t)bucket_char);

        // find the union of all path roots of the bucket
        fmd_tree_t *constraint_values;
        fmd_tree_init(&constraint_values, &prm_fstruct, &prm_fstruct);
        for(int j = 0; j < buckets[i]->size; j++) {
            fmd_graph_ecs_t *path = buckets[i]->data[j];
            fmd_vector_t *incoming_neighbors;
            fmd_table_lookup(graph->incoming_neighbors, (void*)path->head_vid, &incoming_neighbors);
            for(int k = 0; k < incoming_neighbors->size; k++) {
                fmd_tree_insert(constraint_values, incoming_neighbors->data[k], incoming_neighbors->data[k]);
            }
        }
        // flatten keys of the tree
        fmd_vector_t *constraints;
        fmd_vector_init(&constraints, constraint_values->no_items, &prm_fstruct);
        fmd_tree_inorder(constraint_values, constraints, fmd_graph_flatten_constraints);

        // add constraint set to constraint set table, with bucket prefix as key
        fmd_table_insert(constraint_sets, bucket_prefix, constraints);
        fmd_tree_free(constraint_values);

        // call to children recursively
        if(bucket_prefix->size < max_depth) {
            fmd_graph_extract_constraint_sets_helper(buckets[i], bucket_prefix, constraint_sets, alphabet, char2idx,
                                                     idx2char, graph, max_depth);
        }
    }
    // cleanup buckets and paths
    //for(int_t i = 0; i < alphabet->size; i++) {
    //    for(int_t j = 0; j < buckets[i]->size; j++) {
    //        free(buckets[i]->data[j]);
    //    }
    //    fmd_vector_free(buckets[i]);
    //}
    //free(buckets);
}

void fmd_graph_flatten_constraints(void *key, void *value, void *p) {
    fmd_vector_t *constraints = (fmd_vector_t*)p;
    fmd_vector_append(constraints, key);
}
