#ifndef FMD_FMD_DECODER_H
#define FMD_FMD_DECODER_H
#include "fmd_common.h"
#include "fmd_bitstream.h"
#include "fmd_graph.h"

#define FMD_ENCODED_GRAPH_ALPHABET_SIZE_BIT_LENGTH 64
#define FMD_ENCODED_GRAPH_MAX_ENCODING_BIT_LENGTH 8
#define FMD_ENCODED_GRAPH_VID_BIT_LENGTH 64
#define FMD_ENCODED_GRAPH_LABEL_BIT_LENGTH 64
#define FMD_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH 64

typedef struct fmd_encoded_vertex_ {
    fmd_bs_t encoded_string;
    vid_t *outgoing_edges;
    int_t no_outgoing_edges;
    int_t no_encoded_characters;
    vid_t vid;
} fmd_encoded_vertex_t;

typedef struct fmd_encoded_graph_ {
    fmd_encoded_vertex_t *vertices;
    int_t no_vertices;
    int_t no_edges;
    int_t no_total_encoded_characters;
    int_t alphabet_size;
    unsigned char alphabet_occ[256];   // boolean table to see if character is in the alphabet
    unsigned char encoding_table[256]; // each element encodes the log2(alphabet_size) bit encoding value
    unsigned char decoding_table[256]; // each element decodes the encoding value to its alphabet character
} fmd_encoded_graph_t;

void fmd_encoded_graph_init(fmd_encoded_graph_t **encoded_graph, fmd_graph_t *graph);
void fmd_encoded_graph_free(fmd_encoded_graph_t *enc_graph);

/******************************************************************************
 * Serialization and deserialization functions
 *****************************************************************************/

void fmd_encoded_graph_serialize_from_buffer(fmd_encoded_graph_t **encoded_graph, unsigned char *buf, uint64_t buf_size);
/**
 * Write order to buffer:
 *  - FMD_ENCODED_GRAPH_ALPHABET_BIT_LENGTH : alphabet_size
 *  - 256 bytes : alphabet_occ
 *  - 256 bytes : encoding_table
 *  - 256 bytes : decoding_table
 *  - FMD_ENCODED_GRAPH_VID_BIT_LENGTH : no_vertices
 *  - FMD_ENCODED_GRAPH_VID_BIT_LENGTH : no_edges
 *  - FMD_ENCODED_GRAPH_LABEL_BIT_LENGTH : no_total_encoded_characters
 *  - for i = 0 to the number of vertices:
 *      - FMD_ENCODED_GRAPH_VID_BIT_LENGTH : vid
 *      - FMD_ENCODED_GRAPH_LABEL_BIT_LENGTH : no_encoded_characters
 *      - FMD_ENCODED_GRAPH_ADJ_LIST_LEN_BIT_LENGTH : no_outgoing_edges
 *      - for j = 0 to the number of outgoing edges:
 *          - FMD_ENCODED_GRAPH_VID_BIT_LENGTH : edge_vid
 *      - no_encoded_characters * ceil(log2(alphabet_size)) : label payload
 *      - 64-Bit boundary paddding
 */
void fmd_encoded_graph_serialize_to_buffer(fmd_encoded_graph_t *encoded_graph, unsigned char **buf_ret, uint64_t *buf_size_re);

typedef struct fmd_walk_ {
    struct fmd_walk_ *next;
    struct fmd_walk_ *prev;
    void* metadata; // user defined metadata
    vid_t vid; // vid of the substring match
    int_t string_lo;
    int_t string_hi;
    int_t graph_lo;
    int_t graph_hi;
} fmd_walk_node_t;
typedef struct fmd_match_chain_ {
    fmd_walk_node_t *dummy; // dummy->next = head, dummy->prev = tail
    fmd_walk_node_t *head;
    fmd_walk_node_t *tail;
    fmd_encoded_graph_t *graph;
    int_t size;
} fmd_walk_t;
void fmd_walk_init(fmd_walk_t **list, fmd_encoded_graph_t *graph);
void fmd_walk_free(fmd_walk_t *list);
fmd_walk_t *fmd_walk_copy(fmd_walk_t *list);
uint_t fmd_walk_hash(fmd_walk_t *list);
int fmd_walk_comp(fmd_walk_t *l1, fmd_walk_t *l2);

static fmd_fstruct_t fmd_fstruct_walk = {
        (fcomp) fmd_walk_comp,
        (fhash) fmd_walk_hash,
        (ffree) fmd_walk_free,
        (fcopy) fmd_walk_copy,
};

void fmd_walk_append(fmd_walk_t *list, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi);
void fmd_walk_prepend(fmd_walk_t *list, void *metadata, vid_t vid, int_t string_lo, int_t string_hi, int_t graph_lo, int_t graph_hi);

typedef void (*fwalk_extend)(fmd_walk_t *root, fmd_string_t *str, void *metadata, fmd_vector_t **walks);

void fmd_encoded_graph_walk_string(fmd_encoded_graph_t *encoded_graph, fmd_string_t *str, vid_t v, int_t o, void* metadata, fwalk_extend extend_func, fmd_vector_t **walks);
void fmd_encoded_graph_walk_extend_default(fmd_walk_t *root, fmd_string_t *str, void *metadata, fmd_vector_t **tval);


#endif //FMD_FMD_DECODER_H
