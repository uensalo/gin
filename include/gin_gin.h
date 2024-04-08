/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_gin.h is part of gin
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
#ifndef GIN_GIN_GIN_H
#define GIN_GIN_GIN_H
#include "gin_common.h"
#include "gin_graph.h"
#include "gin_interval_merge_tree.h"
#include "gin_oimt.h"
#include "gin_fmi.h"
#include "gin_table.h"
#include "gin_vector.h"
#include "gin_tree.h"
#include "assert.h"

// separator characters
#define GIN_GIN_DEFAULT_c_0 '('
#define GIN_GIN_DEFAULT_c_1 ')'
// permutation encoding characters
#define GIN_GIN_DEFAULT_a_0 ','
#define GIN_GIN_DEFAULT_a_1 '.'
// number of reserved characters
#define GIN_GIN_NO_RESERVED_CHARS 5
// bit field lengths for serialization
#define GIN_GIN_NO_BITS_BIT_LENGTH 64
#define GIN_GIN_SPECIAL_CHAR_BIT_LENGTH 40
#define GIN_GIN_NO_VERTICES_BIT_LENGTH 40
#define GIN_GIN_PERMUTATION_BIT_LENGTH 40
#define GIN_GIN_BWT_TO_VID_BIT_LENGTH 40
#define GIN_GIN_FMI_NO_BITS_BIT_LENGTH 64
#define GIN_GIN_IMT_NO_BITS_BIT_LENGTH 64
#define GIN_GIN_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH 32
#define GIN_GIN_IMT_INTERVAL_BOUNDARY_BIT_LENGTH 40
// bit field lengths for cache
#define GIN_GIN_CACHE_DEPTH_BIT_LENGTH 64
#define GIN_GIN_CACHE_NO_ENTRIES_BIT_LENGTH 64
#define GIN_GIN_CACHE_FMI_SIZE_BIT_LENGTH 64
#define GIN_GIN_CACHE_VALUE_SIZE_BIT_LENGTH 64
#define GIN_GIN_CACHE_FORK_CARDINALITY_BIT_LENGTH 64
#define GIN_GIN_CACHE_FORK_BOUNDARY_BIT_LENGTH 64
#define GIN_GIN_CACHE_FMI_DEFAULT_RANK_RATE 16

#ifdef GIN_SDSL
typedef void* sdsl_csa;
#endif

typedef struct gin_gin_ {
    char_t c_0; // character marking the beginning of a vertex
    char_t c_1; // character marking the end of a vertex
    gin_vector_t *permutation; // permutation to make sa ranges as consecutive as possible
    gin_vector_t *bwt_to_vid; // converts c0 ranks to text ranks, i.e. vids
    int_t *alphabet;
    int_t alphabet_size;
    int_t no_chars;
#ifdef GIN_SDSL
    sdsl_csa *graph_fmi; // fm-index of the graph encoding
#else
    gin_fmi_t *graph_fmi; // fm-index of the graph encoding
#endif
    gin_imt_t *r2r_tree;  // translates sa ranges to sa ranges of incoming nodes
#ifdef GIN_ORACLE
    gin_oimt_t *oracle_r2r;
#endif
} gin_gin_t;

void gin_gin_init(gin_gin_t** gin, gin_graph_t *graph, gin_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate);
gin_vector_t *gin_gin_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords);
void gin_gin_free(gin_gin_t *gin);
int gin_gin_comp(gin_gin_t *f1, gin_gin_t *f2);
void gin_gin_decode(gin_gin_t *gin, gin_graph_t **graph, gin_vector_t **permutation);

// for gathering statistics
typedef struct gin_gin_stats_ {
    int_t no_calls_to_advance_fork;
    int_t no_calls_to_precedence_range;
    int_t no_matching_forks;
    int_t no_partial_forks;
} gin_gin_stats_t;

/******************************************************************************
 * Querying (forks, functions, caching)
 *****************************************************************************/
typedef enum gin_fork_node_type_{
    ROOT = 0,
    MAIN = 1,
    LEAF = 2, // matching leaf
    DEAD = 3, // partial match
    CACH = 4, // cached node, prevents free
} gin_fork_node_type_t;
typedef struct gin_fork_node_{
    int_t sa_lo;
    int_t sa_hi;
    int_t pos;
    gin_fork_node_type_t type;
} gin_fork_node_t;
gin_fork_node_t *gin_fork_node_init(int_t sa_lo, int_t sa_hi,
                                    int_t pos,
                                    gin_fork_node_type_t type);
void gin_fork_node_free(gin_fork_node_t *node);
gin_fork_node_t *gin_fork_node_copy(gin_fork_node_t *node);
uint_t gin_fork_node_hash(gin_fork_node_t *node);
int gin_fork_node_comp(gin_fork_node_t *n1, gin_fork_node_t *n2);
int gin_fork_node_comp_exact(gin_fork_node_t *n1, gin_fork_node_t *n2);

static gin_fstruct_t gin_fstruct_fork_node = {
        (fcomp) gin_fork_node_comp,
        (fhash) gin_fork_node_hash,
        (ffree) gin_fork_node_free,
        (fcopy) gin_fork_node_copy,
};

static gin_fstruct_t gin_fstruct_fork_node_exact = {
        (fcomp) gin_fork_node_comp_exact,
        (fhash) gin_fork_node_hash,
        (ffree) gin_fork_node_free,
        (fcopy) gin_fork_node_copy,
};

typedef struct gin_gin_cache_ { // implements an "FM-table"
    // header entries:
    int_t depth;
    int_t no_entries;
    int_t key_fmi_size_in_bits; // word aligned
    int_t value_buffer_size_in_bits; // word aligned
    word_t *item_offsets;
    word_t *items;
    // payload:
#ifdef GIN_SDSL
    sdsl_csa *key_fmi;
#else
    gin_fmi_t *key_fmi;
#endif
} gin_gin_cache_t;
typedef struct gin_gin_cache_helper_p_ {
    gin_gin_cache_t *cache;
    gin_table_t **cache_tables;
    gin_gin_t *gin;
} gin_gin_cache_helper_p_t;
typedef struct gin_gin_cache_encode_p_ {
    gin_table_t **cache_tables;
    gin_string_t *key_encoding; // the resulting encoding will be written here
    gin_vector_t *values;       // the values will be written here
} gin_gin_cache_encode_p_t;
typedef struct gin_gin_cache_qr_ {
    int_t lo;
    int_t hi;
    int_t pos;
    gin_string_t *pattern;
} gin_gin_cache_qr_t;
void gin_gin_cache_init_step(gin_gin_t *gin, gin_string_t *string, gin_vector_t **cur_forks, gin_vector_t **partial_matches);
void gin_gin_cache_init_helper_trav1(void* key, void* value, void* params); //(*ftrav_kv)(void *key, void *value, void *p);
void gin_gin_cache_init_helper_trav2(void* key, void* value, void* params); //(*ftrav_kv)(void *key, void *value, void *p);
void gin_gin_cache_init(gin_gin_cache_t **cache, gin_gin_t *gin, int_t depth);
bool gin_gin_cache_advance_query(gin_gin_cache_t *fmi, gin_gin_cache_qr_t *qr);
bool gin_gin_cache_query_precedence_range(gin_gin_cache_t *fmi, gin_gin_cache_qr_t *qr, char_t c, int_t *lo, int_t *hi);
void gin_gin_cache_lookup(gin_gin_cache_t *cache, gin_string_t *string, int_t start_pos, int_t max_forks, gin_vector_t **cached_forks);
int_t gin_gin_cache_size(gin_gin_cache_t *cache);

/**
 * Write order to buffer:
 * - GIN_GIN_CACHE_DEPTH_BIT_LENGTH : cache->depth // Number of bits for the entire data structure
 * - GIN_GIN_NO_ENTRIES_BIT_LENGTH: cache->no_entries // Number of <key,value> pairs stored in the cache
 * - for i = 0 to cache->no_entries:
 * -    GIN_GIN_BYTE_OFFSET_BIT_LENGTH: cache->item_offsets[i] // byte offset lookup table into items
 * - cache->key_fmi // FM-index over a special encoding of the keys
 * - cache->items // Item buffer, user interpreted
 *
 * @param cache
 * @param buf_ret
 * @param buf_size_ret
 */
void gin_gin_cache_serialize_to_buffer(gin_gin_cache_t *cache, unsigned char **buf_ret, uint64_t *buf_size_ret);
void gin_gin_cache_serialize_from_buffer(gin_gin_cache_t **cachew, unsigned char *buf, uint64_t buf_size);
void gin_gin_cache_free(gin_gin_cache_t *cache);

bool gin_gin_advance_fork(gin_gin_t *gin, gin_fork_node_t *qr, gin_string_t *pattern);
bool gin_gin_fork_precedence_range(gin_gin_t *gin, gin_fork_node_t *qr, char_t c, int_t *lo, int_t *hi);


// legacy, or for debugging purposes
void gin_gin_query_find_dfs(gin_gin_t *gin, gin_string_t *string, int_t max_forks, gin_vector_t **paths, gin_vector_t **dead_ends, int_t num_threads);
void gin_gin_query_find_dfs_process_fork(gin_gin_t *gin, gin_fork_node_t *fork, int_t max_forks, gin_string_t *pattern, gin_vector_t *exact_matches, gin_vector_t *partial_matches);

void gin_gin_query_find_step(gin_gin_t *gin, gin_string_t *string, int_t max_forks, int_t *t, gin_vector_t **forks, gin_vector_t **partial_matches, gin_gin_stats_t *stats);
void gin_gin_query_find_bootstrapped(gin_gin_t *gin, gin_vector_t *bootstrap, int_t bootstrap_depth, gin_string_t *string, int_t max_forks, gin_vector_t **paths, gin_vector_t **dead_ends, gin_gin_stats_t *stats);
void gin_gin_query_find(gin_gin_t *gin, gin_gin_cache_t *cache, gin_string_t *string, int_t max_forks, gin_vector_t **paths, gin_vector_t **dead_ends, gin_gin_stats_t **stats);
void gin_gin_compact_forks(gin_gin_t *gin, gin_vector_t *forks, gin_vector_t **merged_forks);

/******************************************************************************
 * Result reporting and decoding
 *****************************************************************************/
typedef struct gin_gin_decoded_match_ {
    vid_t vid;
    int_t offset;
} gin_decoded_match_t;

void                 gin_decoded_match_init(gin_decoded_match_t **dec, vid_t vid, int_t offset);
int                  gin_decoded_match_comp(gin_decoded_match_t *dec1, gin_decoded_match_t *dec2);
uint_t               gin_decoded_match_hash(gin_decoded_match_t *dec);
void                 gin_decoded_match_free(gin_decoded_match_t *dec);
gin_decoded_match_t* gin_decoded_match_copy(gin_decoded_match_t *dec);

static gin_fstruct_t gin_fstruct_decoded_match = {
        (fcomp) gin_decoded_match_comp,
        (fhash) gin_decoded_match_hash,
        (ffree) gin_decoded_match_free,
        (fcopy) gin_decoded_match_copy,
};

typedef struct gin_gin_decoder_ {
    gin_gin_t *gin; // fm index for which the decoder will be constructed
    gin_vector_t *vertex_bases;  // vertex base indices
} gin_gin_decoder_t;

void gin_gin_decoder_init(gin_gin_decoder_t **dec, gin_gin_t *gin);
void gin_gin_decoder_free(gin_gin_decoder_t *dec);
/**
 *
 * @param dec GIN decoder object
 * @param sa_lo (inclusive) Start range over the suffix array
 * @param sa_hi (exclusive) End range over the suffix array
 * @param max_matches Number of maximum matches to decode, -1 to decode all
 * @param matches Reference to the output vector containing gin_decoded_match_t objects
 */
void gin_gin_decoder_decode_one(gin_gin_decoder_t *dec, int_t sa_lo, int_t sa_hi, int_t max_matches, gin_vector_t **matches);
/**
 *
 * @param dec GIN decoder object
 * @param matches List of exact matches obtained from a call to one of the gin_query_locate_(...)
 * @param max_matches Number of maximum matches to decode, -1 to decode all
 * @param decoded Reference to the output vector containing gin_vector_t objects containing gin_decoded_match_t objects
 */
void gin_gin_decoder_decode_ends(gin_gin_decoder_t *dec, gin_vector_t *matches, int_t max_matches, gin_vector_t **decoded);

/******************************************************************************
 * Serialization and deserialization functions
 *****************************************************************************/

void gin_gin_serialize_from_buffer(gin_gin_t **gin_ret, unsigned char *buf, uint64_t buf_size);
/**
 * Write order to buffer:
 *  - GIN_GIN_FMI_NO_BITS_BIT_LENGTH : gin_bit_length // Number of bits for the entire data structure
 *  - GIN_GIN_SPECIAL_CHAR_BIT_LENGTH : c_0 // Special character char c_0
 *  - GIN_GIN_SPECIAL_CHAR_BIT_LENGTH : c_1 // Special character char c_1
 *  - GIN_GIN_NO_VERTICES_BIT_LENGTH : no_vertices // number of vertices in the graph of the gin
 *  - for i = 0 to the number of vertices:
 *      - GIN_GIN_PERMUTATION_BIT_LENGTH permutation[i] // Permutation entry
 *  - for i = 0 to the number of vertices:
 *      - GIN_GIN_BWT_TO_VID_BIT_LENGTH bwt_to_vid[i] // Rank translation entry
 *  - GIN_GIN_FMI_NO_BITS_BIT_LENGTH : fmi_bit_length // bit-length of the FMI bitstream
 *  - fmi_bit_length : fmi // the buffer for the fm-index, word aligned
 *  - GIN_GIN_IMT_NO_BITS_BIT_LENGTH : imt_bit_length // bit_length of the IMT buffer
 *  - imt_bit_length: imt // imt buffer
 * @param gin
 * @param buf
 * @param buf_size
 */
void gin_gin_serialize_to_buffer(gin_gin_t *gin, unsigned char **buf_ret, uint64_t *buf_size_re);
void gin_gin_serialize_to_buffer_imt_helper(gin_imt_node_t *node, gin_bs_t *bs, uint_t *widx);


#endif //GIN_GIN_GIN_H
