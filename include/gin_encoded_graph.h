#ifndef GIN_GIN_DECODER_H
#define GIN_GIN_DECODER_H
#include "gin_common.h"
#include "gin_bitstream.h"
#include "gin_graph.h"

#define GIN_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH 64
#define GIN_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH 8
#define GIN_ENCODED_GRAPH_VID_BIT_LENGTH 64
#define GIN_ENCODED_GRAPH_LABEL_BIT_LENGTH 64
#define GIN_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH 64

typedef struct gin_encoded_vertex_ {
    gin_bs_t encoded_string;
    vid_t *outgoing_edges;
    int_t no_outgoing_edges;
    int_t no_encoded_characters;
    vid_t vid;
} gin_encoded_vertex_t;

typedef struct gin_encoded_graph_ {
    gin_encoded_vertex_t *vertices;
    int_t no_vertices;
    int_t no_edges;
    int_t no_total_encoded_characters;
    int_t alphabet_size;
    unsigned char alphabet_occ[256];   // boolean table to see if character is in the alphabet
    unsigned char encoding_table[256]; // each element encodes the log2(alphabet_size) bit encoding value
    unsigned char decoding_table[256]; // each element decodes the encoding value to its alphabet character
} gin_encoded_graph_t;

void gin_encoded_graph_init(gin_encoded_graph_t **encoded_graph, gin_graph_t *graph);
void gin_encoded_graph_free(gin_encoded_graph_t *enc_graph);

/******************************************************************************
 * Serialization and deserialization functions
 *****************************************************************************/

void gin_encoded_graph_serialize_from_buffer(gin_encoded_graph_t **encoded_graph, unsigned char *buf, uint64_t buf_size);
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
void gin_encoded_graph_serialize_to_buffer(gin_encoded_graph_t *encoded_graph, unsigned char **buf_ret, uint64_t *buf_size_re);

typedef struct gin_walk_ {
    struct gin_walk_ *next;
    struct gin_walk_ *prev;
    void* metadata; // user defined metadata
    vid_t vid; // vid of the substring match
    int_t string_lo;
    int_t string_hi;
    int_t graph_lo;
    int_t graph_hi;
} gin_walk_node_t;
typedef struct gin_match_chain_ {
    gin_walk_node_t *dummy; // dummy->next = head, dummy->prev = tail
    gin_walk_node_t *head;
    gin_walk_node_t *tail;
    gin_encoded_graph_t *graph;
    int_t size;
} gin_walk_t;
void gin_walk_init(gin_walk_t **list, gin_encoded_graph_t *graph);
void gin_walk_free(gin_walk_t *list);
gin_walk_t *gin_walk_copy(gin_walk_t *list);
uint_t gin_walk_hash(gin_walk_t *list);
int gin_walk_comp(gin_walk_t *l1, gin_walk_t *l2);

static gin_fstruct_t gin_fstruct_walk = {
        (fcomp) gin_walk_comp,
        (fhash) gin_walk_hash,
        (ffree) gin_walk_free,
        (fcopy) gin_walk_copy,
};

void gin_walk_append(gin_walk_t *list, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi);
void gin_walk_prepend(gin_walk_t *list, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi);

typedef void (*fwalk_extend)(gin_walk_t *root, gin_string_t *str, void *metadata, gin_vector_t **walks);

void gin_encoded_graph_walk_string(gin_encoded_graph_t *encoded_graph, gin_string_t *str, vid_t v, int_t o, void* metadata, fwalk_extend extend_func, gin_vector_t **walks);
void gin_encoded_graph_walk_extend_default(gin_walk_t *root, gin_string_t *str, void *metadata, gin_vector_t **tval);


#endif //GIN_GIN_DECODER_H
