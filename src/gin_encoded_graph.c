#include "gin_encoded_graph.h"

void gin_encoded_graph_init(gin_encoded_graph_t **encoded_graph, gin_graph_t *graph) {
    gin_encoded_graph_t *eg = calloc(1, sizeof(gin_encoded_graph_t));
    if(!eg) return;
    eg->no_vertices = graph->vertex_list->size;
    eg->vertices = calloc(eg->no_vertices, sizeof(gin_encoded_vertex_t));
    // gather the alphabet and construct adjacency lists
    for(int_t i = 0; i < graph->vertex_list->size; i++) {
        gin_vertex_t *vertex = graph->vertex_list->data[i];
        int_t vid = vertex->id;
        // create the adjacency list for the encoded graph
        eg->vertices[vid].vid = vid;
        eg->vertices[vid].no_encoded_characters = vertex->label->size;
        void *_;
        gin_table_lookup(graph->outgoing_neighbors, (void*)vid, &_);
        gin_vector_t *outgoing_edges = (gin_vector_t*)_;
        eg->vertices[vid].no_outgoing_edges = outgoing_edges->size;
        eg->vertices[vid].outgoing_edges = calloc(outgoing_edges->size, sizeof(vid_t));
        for(int_t j = 0; j < outgoing_edges->size; j++) {
            eg->vertices[vid].outgoing_edges[j] = (vid_t)outgoing_edges->data[j];
        }
        for(int_t j = 0; j < vertex->label->size; j++) {
            eg->alphabet_occ[vertex->label->seq[j]] = 1; // mark as existing
        }
    }
    // construct the encoding and decoding tables according to alphabet occupancy
    int_t k = 0;
    for(int_t i = 0; i < 256; i++) {
        if(eg->alphabet_occ[i]) {
            eg->encoding_table[i] = k;
            eg->decoding_table[k] = (char)i;
            ++eg->alphabet_size;
            ++k;
        }
    }
    // encode vertex labels
    int_t bits_per_char = gin_ceil_log2(eg->alphabet_size);
    for(int_t i = 0; i < graph->vertex_list->size; i++) {
        gin_vertex_t *vertex = graph->vertex_list->data[i];
        vid_t vid = vertex->id;
        gin_bs_t *bs;
        int_t no_words_in_encoding = 1 + (vertex->label->size * bits_per_char - 1) / WORD_NUM_BITS;
        gin_bs_init_reserve(&bs, no_words_in_encoding);
        int_t idx = 0;
        for(int_t j = 0; j < vertex->label->size; j++) {
            char_t c = vertex->label->seq[j];
            gin_bs_write_word(bs, idx, (word_t)eg->encoding_table[c], bits_per_char);
            idx += bits_per_char;
        }
        eg->vertices[vid].encoded_string = *bs;
        gin_bs_free_disown(bs);
    }

    // gather batch statistics
    for(int_t i = 0; i < eg->no_vertices; i++) {
        eg->no_total_encoded_characters += eg->vertices[i].no_encoded_characters;
        eg->no_edges += eg->vertices[i].no_outgoing_edges;
    }
    *encoded_graph = eg;
}

void gin_encoded_graph_free(gin_encoded_graph_t *enc_graph) {
    if(!enc_graph) return;
    for(int_t i = 0; i < enc_graph->no_vertices; i++) {
        free(enc_graph->vertices[i].outgoing_edges);
        gin_bs_free_no_alloc(&enc_graph->vertices[i].encoded_string);
    }
    free(enc_graph->vertices);
    free(enc_graph);
}

/******************************************************************************
 * Serialization and deserialization functions
 *****************************************************************************/

void gin_encoded_graph_serialize_from_buffer(gin_encoded_graph_t **encoded_graph, unsigned char *buf, uint64_t buf_size) {
    gin_encoded_graph_t *eg = calloc(1, sizeof(gin_encoded_graph_t));
    gin_bs_t *bs;
    gin_bs_init_from_buffer(buf, buf_size, &bs);
    word_t ridx = 0;

    word_t alphabet_size = 0;
    gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH, &alphabet_size);
    ridx += GIN_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH;
    eg->alphabet_size = (int_t)alphabet_size;

    int_t bits_per_char = gin_ceil_log2(alphabet_size);

    for(int_t i = 0; i < 256; i++) {
        word_t occ = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH, &occ);
        ridx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
        eg->alphabet_occ[i] = occ;
    }
    // write encoding_table
    for(int_t i = 0; i < 256; i++) {
        word_t enc = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH, &enc);
        ridx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
        eg->encoding_table[i] = enc;
    }
    // write decoding table
    for(int_t i = 0; i < 256; i++) {
        word_t dec = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH, &dec);
        ridx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
        eg->decoding_table[i] = dec;
    }

    word_t no_vertices = 0;
    gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_VID_BIT_LENGTH, &no_vertices);
    ridx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;
    eg->no_vertices = (int_t)no_vertices;

    eg->vertices = calloc(eg->no_vertices, sizeof(gin_encoded_vertex_t));

    word_t no_edges = 0;
    gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_VID_BIT_LENGTH, &no_edges);
    ridx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;
    eg->no_edges = (int_t)no_edges;

    word_t no_total_encoded_characters = 0;
    gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH, &no_total_encoded_characters);
    ridx += GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH;
    eg->no_total_encoded_characters = (int_t)no_total_encoded_characters;

    for(int_t i = 0; i < no_vertices; i++) {
        word_t rvid = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_VID_BIT_LENGTH, &rvid);
        ridx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;
        vid_t vid = (vid_t)rvid;

        word_t no_encoded_characters = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH, &no_encoded_characters);
        ridx += GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH;
        eg->vertices[vid].no_encoded_characters = (int_t)no_encoded_characters;

        word_t no_outgoing_edges = 0;
        gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH, &no_outgoing_edges);
        ridx += GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH;
        eg->vertices[vid].no_outgoing_edges = (int_t)no_outgoing_edges;

        // alloc adjacency list
        eg->vertices[vid].outgoing_edges = calloc(no_outgoing_edges,sizeof(vid_t));

        for(int_t j = 0 ; j < eg->vertices[vid].no_outgoing_edges; j++) {
            word_t edge_vid = 0;
            gin_bs_read_word(bs, ridx, GIN_ENCODED_GRAPH_VID_BIT_LENGTH, &edge_vid);
            ridx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;
            eg->vertices[vid].outgoing_edges[j] = (int_t)edge_vid;
        }
        int_t no_encoded_characters_bits = bits_per_char * (int_t)no_encoded_characters;
        int_t no_encoded_words = 1 + ((no_encoded_characters_bits - 1) >> WORD_LOG_BITS);

        gin_bs_init_reserve_no_alloc(&eg->vertices[vid].encoded_string, no_encoded_words);

        word_t widx = 0;
        for(int_t j = 0; j < no_encoded_words; j++) {
            word_t encoded_bits = 0;
            gin_bs_read_word(bs, ridx, WORD_NUM_BITS, &encoded_bits);
            ridx += WORD_NUM_BITS;
            gin_bs_write_word(&eg->vertices[vid].encoded_string, widx, encoded_bits, WORD_NUM_BITS);
            widx += WORD_NUM_BITS;
        }
        ridx = (1 + ((ridx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    }
    gin_bs_free(bs);
    *encoded_graph = eg;
}
/**
 * Write order to buffer:
 *  - GIN_ENCODED_GRAPH_ALPHABET_BIT_LENGTH : alphabet_size
 *  - 256 bytes : alphabet_occ
 *  - 256 bytes : encoding_table
 *  - 256 bytes : decoding_table
 *  - GIN_ENCODED_GRAPH_VID_BIT_LENGTH : no_vertices
 *  - GIN_ENCODED_GRAPH_VID_BIT_LENGTH : no_edges
 *  - GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH : no_total_encoded_characters
 *  - for i = 0 to the number of vertices:
 *      - GIN_ENCODED_GRAPH_VID_BIT_LENGTH : vid
 *      - GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH : no_encoded_characters
 *      - GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH : no_outgoing_edges
 *      - for j = 0 to the number of outgoing edges:
 *          - GIN_ENCODED_GRAPH_VID_BIT_LENGTH : edge_vid
 *      - no_encoded_characters * ceil(log2(alphabet_size)) : label payload
 *      - 64-Bit boundary paddding
 */
void gin_encoded_graph_serialize_to_buffer(gin_encoded_graph_t *encoded_graph, unsigned char **buf_ret, uint64_t *buf_size_ret) {
    gin_bs_t *bits;
    gin_bs_init(&bits);
    word_t widx = 0;
    int_t bits_per_char = gin_ceil_log2(encoded_graph->alphabet_size);

    gin_bs_write_word(bits, widx, encoded_graph->alphabet_size, GIN_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH);
    widx +=  GIN_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH;
    // write alphabet occ
    for(int_t i = 0; i < 256; i++) {
        word_t occ = (word_t)encoded_graph->alphabet_occ[i];
        gin_bs_write_word(bits, widx, occ, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
    }
    // write encoding_table
    for(int_t i = 0; i < 256; i++) {
        word_t enc = (word_t)encoded_graph->encoding_table[i];
        gin_bs_write_word(bits, widx, enc, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
    }
    // write decoding table
    for(int_t i = 0; i < 256; i++) {
        word_t dec = (word_t)encoded_graph->decoding_table[i];
        gin_bs_write_word(bits, widx, dec, GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH;
    }

    gin_bs_write_word(bits, widx, encoded_graph->no_vertices, GIN_ENCODED_GRAPH_VID_BIT_LENGTH);
    widx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;

    gin_bs_write_word(bits, widx, encoded_graph->no_edges, GIN_ENCODED_GRAPH_VID_BIT_LENGTH);
    widx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;

    gin_bs_write_word(bits, widx, encoded_graph->no_total_encoded_characters, GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH);
    widx += GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH;

    for(int_t i = 0; i < encoded_graph->no_vertices; i++) {
        gin_bs_write_word(bits, widx, encoded_graph->vertices[i].vid, GIN_ENCODED_GRAPH_VID_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;

        gin_bs_write_word(bits, widx, encoded_graph->vertices[i].no_encoded_characters, GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH;

        gin_bs_write_word(bits, widx, encoded_graph->vertices[i].no_outgoing_edges, GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH);
        widx += GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH;

        for(int_t j = 0; j < encoded_graph->vertices[i].no_outgoing_edges; j++) {
            gin_bs_write_word(bits, widx, encoded_graph->vertices[i].outgoing_edges[j], GIN_ENCODED_GRAPH_VID_BIT_LENGTH);
            widx += GIN_ENCODED_GRAPH_VID_BIT_LENGTH;
        }
        for(int_t j = 0; j < encoded_graph->vertices[i].encoded_string.cap_in_words; j++) {
            gin_bs_write_word(bits, widx, encoded_graph->vertices[i].encoded_string.words[j], WORD_NUM_BITS);
            widx += WORD_NUM_BITS;
        }
        // align to word boundary
        widx = (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    }
    word_t* buf;
    uint64_t no_words;
    gin_bs_fit(bits, widx);
    gin_bs_detach(bits, &buf, &no_words);
    gin_bs_free(bits);
    *buf_ret = (unsigned char*)buf;
    *buf_size_ret = no_words * sizeof(word_t);
}


void gin_walk_init(gin_walk_t **walk, gin_encoded_graph_t *graph) {
    *walk = (gin_walk_t*)calloc(1, sizeof(gin_walk_t));
    (*walk)->dummy = (gin_walk_node_t*)calloc(1, sizeof(gin_walk_node_t));
    (*walk)->dummy->next = (*walk)->dummy;
    (*walk)->dummy->prev = (*walk)->dummy;
    (*walk)->graph = graph;
    (*walk)->head = NULL;
    (*walk)->tail = NULL;
    (*walk)->size = 0;
}

void gin_walk_free(gin_walk_t *walk) {
    gin_walk_node_t *current = walk->dummy->next;
    while(current != walk->dummy) {
        gin_walk_node_t *next = current->next;
        free(current);
        current = next;
    }
    free(walk->dummy);
    free(walk);
}

gin_walk_t *gin_walk_copy(gin_walk_t *walk) {
    gin_walk_t *copy;
    gin_walk_init(&copy, walk->graph);

    gin_walk_node_t *current = walk->dummy->next;
    while(current != walk->dummy) {
        gin_walk_append(copy, NULL, current->vid, current->string_lo, current->string_hi, current->graph_lo, current->graph_hi);
        current = current->next;
    }
    return copy;
}

uint_t gin_walk_hash(gin_walk_t *walk) {
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;
    gin_walk_node_t *current = walk->dummy->next;
    while(current != walk->dummy) {
        uint_t hash_item = prm_hash_f((void*)current->vid) ^
                           prm_hash_f((void*)current->string_lo) ^
                           prm_hash_f((void*)current->string_hi) ^
                           prm_hash_f((void*)current->graph_lo) ^
                           prm_hash_f((void*)current->graph_hi);
        hash ^= hash_item;
        hash *= prime;
        current = current->next;
    }
    return hash;
}

int gin_walk_comp(gin_walk_t *l1, gin_walk_t *l2) {
    gin_walk_node_t *n1 = l1->dummy->next;
    gin_walk_node_t *n2 = l2->dummy->next;
    while(n1 != l1->dummy && n2 != l2->dummy) {
        int cmp = n1->vid == n2->vid &&
                  n1->string_lo == n2->string_lo &&
                  n1->string_hi == n2->string_hi &&
                  n1->graph_lo == n2->graph_lo &&
                  n1->graph_hi == n2->graph_hi;
        if(cmp != 0) {
            return cmp;
        }
        n1 = n1->next;
        n2 = n2->next;
    }
    return (n1 != l1->dummy) - (n2 != l2->dummy);
}

void gin_walk_append(gin_walk_t *walk, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi) {
    gin_walk_node_t *new_node = (gin_walk_node_t*)calloc(1,sizeof(gin_walk_node_t));
    new_node->metadata = metadata;
    new_node->vid = vid;
    new_node->string_lo = string_lo;
    new_node->string_hi = string_hi;
    new_node->graph_lo = graph_lo;
    new_node->graph_hi = graph_hi;
    if (walk->size == 0) {
        walk->head = new_node;
        walk->tail = new_node;
        new_node->prev = walk->dummy;
        new_node->next = walk->dummy;
        walk->dummy->prev = new_node;
        walk->dummy->next = new_node;
    } else {
        new_node->prev = walk->tail;
        new_node->next = walk->dummy;
        walk->tail->next = new_node;
        walk->dummy->prev = new_node;
        walk->tail = new_node;
    }
    walk->size++;
}

void gin_walk_prepend(gin_walk_t *walk, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi) {
    gin_walk_node_t *new_node = (gin_walk_node_t*)calloc(1,sizeof(gin_walk_node_t));
    new_node->metadata = metadata;
    new_node->vid = vid;
    new_node->string_lo = string_lo;
    new_node->string_hi = string_hi;
    new_node->graph_lo = graph_lo;
    new_node->graph_hi = graph_hi;
    if (walk->size == 0) {
        walk->head = new_node;
        walk->tail = new_node;
        new_node->prev = walk->dummy;
        new_node->next = walk->dummy;
        walk->dummy->prev = new_node;
        walk->dummy->next = new_node;
    } else {
        new_node->next = walk->head;
        new_node->prev = walk->dummy;
        walk->head->prev = new_node;
        walk->dummy->next = new_node;
        walk->head = new_node;
    }
    walk->size++;
}

void gin_encoded_graph_walk_string(gin_encoded_graph_t *encoded_graph, gin_string_t *str, vid_t v, int_t o, void* metadata, fwalk_extend extend_func, gin_vector_t **walks_ret) {
    gin_vector_t *walks;
    gin_vector_init(&walks, GIN_VECTOR_INIT_SIZE, &gin_fstruct_walk);
    // init the root as the start of the walk
    gin_walk_t *root;
    gin_walk_init(&root, encoded_graph);
    gin_walk_append(root, NULL, v, 0, 0, o, o);
    // extend
    extend_func(root, str, metadata, &walks);
    *walks_ret = walks;
}

void gin_encoded_graph_walk_extend_default(gin_walk_t *root, gin_string_t *str, void *metadata, gin_vector_t **walks) {
    gin_bs_t *enc_str = (gin_bs_t*)metadata;
    int_t bits_per_char = gin_ceil_log2(root->graph->alphabet_size);
    // get the last inserted node
    gin_walk_node_t *last = root->tail;
    gin_encoded_vertex_t *v = &root->graph->vertices[last->vid];
    int_t no_chars_to_match = MIN2(str->size - last->string_hi, v->no_encoded_characters - last->graph_hi);

    // check if the string matches
    int_t encoding_length = (no_chars_to_match * bits_per_char);
    int_t no_words = encoding_length >> WORD_LOG_BITS;
    int_t no_slack_bits = encoding_length & WORD_LOG_MASK;
    word_t sidx = bits_per_char * last->string_hi;
    word_t gidx = bits_per_char * last->graph_hi;

    for(int_t i = 0; i < no_words; i++) {
        word_t sword = 0;
        word_t gword = 0;
        gin_bs_read_word(enc_str, sidx, WORD_NUM_BITS, &sword);
        gin_bs_read_word(&v->encoded_string, gidx, WORD_NUM_BITS, &gword);
        if(sword != gword)  {
            gin_walk_free(root);
            return;
        }
        sidx += WORD_NUM_BITS;
        gidx += WORD_NUM_BITS;
    }
    if(no_slack_bits) {
        word_t sword = 0;
        word_t gword = 0;
        gin_bs_read_word(enc_str, sidx, no_slack_bits, &sword);
        gin_bs_read_word(&v->encoded_string, gidx, no_slack_bits, &gword);
        if(sword != gword) {
            gin_walk_free(root);
            return;
        }
    }
    // mark match on the last node
    last->graph_hi += no_chars_to_match;
    last->string_hi += no_chars_to_match;

    bool string_exhausted = last->string_hi == str->size;
    bool vertex_exhausted = last->graph_hi == v->no_encoded_characters;

    // extend the walk if string is not exhausted
    if(vertex_exhausted && !string_exhausted) {
        if(v->outgoing_edges) {
            for (int_t i = 1; i < v->no_outgoing_edges; i++) {
                gin_walk_t *copy = gin_walk_copy(root);
                gin_walk_append(copy, NULL, v->outgoing_edges[i], last->string_hi, last->string_hi, 0, 0);
                gin_encoded_graph_walk_extend_default(copy, str, metadata, walks);
            }
            gin_walk_append(root, NULL, v->outgoing_edges[0], last->string_hi, last->string_hi, 0, 0);
            gin_encoded_graph_walk_extend_default(root, str, metadata, walks);
        } else {
            gin_walk_free(root);
        }
    }
    else if(string_exhausted) {
        gin_vector_append(*walks, root);
    } else {
        //gin_walk_free(root);
    }
}