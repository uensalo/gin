/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_fmd.h is part of fmd
 *
 * fmd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * fmd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FMD_FMD_FMD_H
#define FMD_FMD_FMD_H
#include "fmd_common.h"
#include "fmd_graph.h"
#include "fmd_interval_merge_tree.h"
#include "fmd_table.h"
#include "fmd_vector.h"
#include "fmd_tree.h"
#include "fmd_bitstream.h"
#include "fmd_fmi.h"
#include "assert.h"

// separator characters
#define FMD_FMD_DEFAULT_c_0 '('
#define FMD_FMD_DEFAULT_c_1 ')'
// permutation encoding characters
#define FMD_FMD_DEFAULT_a_0 ','
#define FMD_FMD_DEFAULT_a_1 '.'
// number of reserved characters
#define FMD_FMD_NO_RESERVED_CHARS 5
// bit field lengths for serialization
#define FMD_FMD_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_SPECIAL_CHAR_BIT_LENGTH 40
#define FMD_FMD_NO_VERTICES_BIT_LENGTH 40
#define FMD_FMD_PERMUTATION_BIT_LENGTH 40
#define FMD_FMD_BWT_TO_VID_BIT_LENGTH 40
#define FMD_FMD_FMI_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_IMT_NO_BITS_BIT_LENGTH 64
#define FMD_FMD_IMT_INTERVAL_LIST_LENGTH_BIT_LENGTH 32
#define FMD_FMD_IMT_INTERVAL_BOUNDARY_BIT_LENGTH 40
// bit field lengths for cache
#define FMD_FMD_CACHE_DEPTH_BIT_LENGTH 64
#define FMD_FMD_CACHE_NO_ENTRIES_BIT_LENGTH 64
#define FMD_FMD_CACHE_FMI_SIZE_BIT_LENGTH 64
#define FMD_FMD_CACHE_VALUE_SIZE_BIT_LENGTH 64
#define FMD_FMD_CACHE_FORK_CARDINALITY_BIT_LENGTH 64
#define FMD_FMD_CACHE_FORK_BOUNDARY_BIT_LENGTH 64
#define FMD_FMD_CACHE_FMI_DEFAULT_RANK_RATE 16

typedef void* sdsl_csa;

typedef struct fmd_fmd_ {
    char_t c_0; // character marking the beginning of a vertex
    char_t c_1; // character marking the end of a vertex
    fmd_vector_t *permutation; // permutation to make sa ranges as consecutive as possible
    fmd_vector_t *bwt_to_vid; // converts c0 ranks to text ranks, i.e. vids
    sdsl_csa *graph_fmi; // fm-index of the graph encoding
    fmd_imt_t *r2r_tree;  // translates sa ranges to sa ranges of incoming nodes
} fmd_fmd_t;

void fmd_fmd_init(fmd_fmd_t** fmd, fmd_graph_t *graph, fmd_vector_t *permutation, char_t c_0, char_t c_1, int_t rank_sample_rate, int_t isa_sample_rate);
fmd_vector_t *fmd_fmd_init_pcodes_fixed_binary_helper(char_t a_0, char_t a_1, int_t no_codewords);
void fmd_fmd_free(fmd_fmd_t *fmd);
int fmd_fmd_comp(fmd_fmd_t *f1, fmd_fmd_t *f2);

// for gathering statistics
typedef struct fmd_fmd_stats_ {
    int_t no_calls_to_advance_fork;
    int_t no_calls_to_precedence_range;
    int_t no_matching_forks;
    int_t no_partial_forks;
} fmd_fmd_stats_t;

/******************************************************************************
 * Querying (forks, functions, caching)
 *****************************************************************************/
typedef enum fmd_fork_node_type_{
    ROOT = 0,
    MAIN = 1,
    LEAF = 2, // matching leaf
    DEAD = 3, // partial match
    CACH = 4, // cached node, prevents free
} fmd_fork_node_type_t;
typedef struct fmd_fork_node_{
    int_t sa_lo;
    int_t sa_hi;
    int_t pos;
    fmd_fork_node_type_t type;
} fmd_fork_node_t;
fmd_fork_node_t *fmd_fork_node_init(int_t sa_lo, int_t sa_hi,
                                    int_t pos,
                                    fmd_fork_node_type_t type);
void fmd_fork_node_free(fmd_fork_node_t *node);
fmd_fork_node_t *fmd_fork_node_copy(fmd_fork_node_t *node);
uint_t fmd_fork_node_hash(fmd_fork_node_t *node);
int fmd_fork_node_comp(fmd_fork_node_t *n1, fmd_fork_node_t *n2);
int fmd_fork_node_comp_exact(fmd_fork_node_t *n1, fmd_fork_node_t *n2);

static fmd_fstruct_t fmd_fstruct_fork_node = {
        (fcomp) fmd_fork_node_comp,
        (fhash) fmd_fork_node_hash,
        (ffree) fmd_fork_node_free,
        (fcopy) fmd_fork_node_copy,
};

static fmd_fstruct_t fmd_fstruct_fork_node_exact = {
        (fcomp) fmd_fork_node_comp_exact,
        (fhash) fmd_fork_node_hash,
        (ffree) fmd_fork_node_free,
        (fcopy) fmd_fork_node_copy,
};

typedef struct fmd_fmd_cache_ { // implements an "FM-table"
    // header entries:
    int_t depth;
    int_t no_entries;
    int_t key_fmi_size_in_bits; // word aligned
    int_t value_buffer_size_in_bits; // word aligned
    word_t *item_offsets;
    // payload:
    fmd_fmi_t *key_fmi;
    word_t *items;
    // not stored, derived
    unsigned char *disk_buffer;
} fmd_fmd_cache_t;
typedef struct fmd_fmd_cache_helper_p_ {
    fmd_fmd_cache_t *cache;
    fmd_table_t **cache_tables;
    fmd_fmd_t *fmd;
} fmd_fmd_cache_helper_p_t;
typedef struct fmd_fmd_cache_encode_p_ {
    fmd_table_t **cache_tables;
    fmd_string_t *key_encoding; // the resulting encoding will be written here
    fmd_vector_t *values;       // the values will be written here
} fmd_fmd_cache_encode_p_t;
void fmd_fmd_cache_init_step(fmd_fmd_t *fmd, fmd_string_t *string, fmd_vector_t **cur_forks, fmd_vector_t **partial_matches);
void fmd_fmd_cache_init_helper_trav1(void* key, void* value, void* params); //(*ftrav_kv)(void *key, void *value, void *p);
void fmd_fmd_cache_init_helper_trav2(void* key, void* value, void* params); //(*ftrav_kv)(void *key, void *value, void *p);
void fmd_fmd_cache_init(fmd_fmd_cache_t **cache, fmd_fmd_t *fmd, int_t depth);
void fmd_fmd_cache_lookup(fmd_fmd_cache_t *cache, fmd_string_t *string, int_t start_pos, int_t max_forks, fmd_vector_t **cached_forks);
int_t fmd_fmd_cache_size(fmd_fmd_cache_t *cache);

/**
 * Write order to buffer:
 * - FMD_FMD_CACHE_DEPTH_BIT_LENGTH : cache->depth // Number of bits for the entire data structure
 * - FMD_FMD_NO_ENTRIES_BIT_LENGTH: cache->no_entries // Number of <key,value> pairs stored in the cache
 * - for i = 0 to cache->no_entries:
 * -    FMD_FMD_BYTE_OFFSET_BIT_LENGTH: cache->item_offsets[i] // byte offset lookup table into items
 * - cache->key_fmi // FM-index over a special encoding of the keys
 * - cache->items // Item buffer, user interpreted
 *
 * @param cache
 * @param buf_ret
 * @param buf_size_ret
 */
void fmd_fmd_cache_serialize_to_buffer(fmd_fmd_cache_t *cache, unsigned char **buf_ret, uint64_t *buf_size_ret);
void fmd_fmd_cache_serialize_from_buffer(fmd_fmd_cache_t **cachew, unsigned char *buf, uint64_t buf_size);
void fmd_fmd_cache_free(fmd_fmd_cache_t *cache);

bool fmd_fmd_advance_fork(fmd_fmd_t *fmd, fmd_fork_node_t *qr, fmd_string_t *pattern);
bool fmd_fmd_fork_precedence_range(fmd_fmd_t *fmd, fmd_fork_node_t *qr, char_t c, int_t *lo, int_t *hi);

// legacy, or for debugging purposes
void fmd_fmd_query_find_dfs(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends, int_t num_threads);
void fmd_fmd_query_find_dfs_process_fork(fmd_fmd_t *fmd, fmd_fork_node_t *fork, int_t max_forks, fmd_string_t *pattern, fmd_vector_t *exact_matches, fmd_vector_t *partial_matches);

void fmd_fmd_query_find_step(fmd_fmd_t *fmd, fmd_string_t *string, int_t max_forks, int_t *t, fmd_vector_t **forks, fmd_vector_t **partial_matches, fmd_fmd_stats_t *stats);
void fmd_fmd_query_find_bootstrapped(fmd_fmd_t *fmd, fmd_vector_t *bootstrap, int_t bootstrap_depth, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends, fmd_fmd_stats_t *stats);
void fmd_fmd_query_find(fmd_fmd_t *fmd, fmd_fmd_cache_t *cache, fmd_string_t *string, int_t max_forks, fmd_vector_t **paths, fmd_vector_t **dead_ends, fmd_fmd_stats_t **stats);
void fmd_fmd_compact_forks(fmd_fmd_t *fmd, fmd_vector_t *forks, fmd_vector_t **merged_forks);

/******************************************************************************
 * Result reporting and decoding
 *****************************************************************************/
typedef struct fmd_fmd_decoded_match_ {
    vid_t vid;
    int_t offset;
} fmd_decoded_match_t;

void                 fmd_decoded_match_init(fmd_decoded_match_t **dec, vid_t vid, int_t offset);
int                  fmd_decoded_match_comp(fmd_decoded_match_t *dec1, fmd_decoded_match_t *dec2);
uint_t               fmd_decoded_match_hash(fmd_decoded_match_t *dec);
void                 fmd_decoded_match_free(fmd_decoded_match_t *dec);
fmd_decoded_match_t* fmd_decoded_match_copy(fmd_decoded_match_t *dec);

static fmd_fstruct_t fmd_fstruct_decoded_match = {
        (fcomp) fmd_decoded_match_comp,
        (fhash) fmd_decoded_match_hash,
        (ffree) fmd_decoded_match_free,
        (fcopy) fmd_decoded_match_copy,
};

typedef struct fmd_fmd_decoder_ {
    fmd_fmd_t *fmd; // fm index for which the decoder will be constructed
    fmd_vector_t *vertex_bases;  // vertex base indices
} fmd_fmd_decoder_t;

void fmd_fmd_decoder_init(fmd_fmd_decoder_t **dec, fmd_fmd_t *fmd);
void fmd_fmd_decoder_free(fmd_fmd_decoder_t *dec);
/**
 *
 * @param dec FMD decoder object
 * @param sa_lo (inclusive) Start range over the suffix array
 * @param sa_hi (exclusive) End range over the suffix array
 * @param max_matches Number of maximum matches to decode, -1 to decode all
 * @param matches Reference to the output vector containing fmd_decoded_match_t objects
 */
void fmd_fmd_decoder_decode_one(fmd_fmd_decoder_t *dec, int_t sa_lo, int_t sa_hi, int_t max_matches, fmd_vector_t **matches);
/**
 *
 * @param dec FMD decoder object
 * @param matches List of exact matches obtained from a call to one of the fmd_query_locate_(...)
 * @param max_matches Number of maximum matches to decode, -1 to decode all
 * @param decoded Reference to the output vector containing fmd_vector_t objects containing fmd_decoded_match_t objects
 */
void fmd_fmd_decoder_decode_ends(fmd_fmd_decoder_t *dec, fmd_vector_t *matches, int_t max_matches, fmd_vector_t **decoded);

/******************************************************************************
 * Serialization and deserialization functions
 *****************************************************************************/

void fmd_fmd_serialize_from_buffer(fmd_fmd_t **fmd_ret, unsigned char *buf, uint64_t buf_size);
/**
 * Write order to buffer:
 *  - FMD_FMD_FMI_NO_BITS_BIT_LENGTH : fmd_bit_length // Number of bits for the entire data structure
 *  - FMD_FMD_SPECIAL_CHAR_BIT_LENGTH : c_0 // Special character char c_0
 *  - FMD_FMD_SPECIAL_CHAR_BIT_LENGTH : c_1 // Special character char c_1
 *  - FMD_FMD_NO_VERTICES_BIT_LENGTH : no_vertices // number of vertices in the graph of the fmd
 *  - for i = 0 to the number of vertices:
 *      - FMD_FMD_PERMUTATION_BIT_LENGTH permutation[i] // Permutation entry
 *  - for i = 0 to the number of vertices:
 *      - FMD_FMD_BWT_TO_VID_BIT_LENGTH bwt_to_vid[i] // Rank translation entry
 *  - FMD_FMD_FMI_NO_BITS_BIT_LENGTH : fmi_bit_length // bit-length of the FMI bitstream
 *  - fmi_bit_length : fmi // the buffer for the fm-index, word aligned
 *  - FMD_FMD_IMT_NO_BITS_BIT_LENGTH : imt_bit_length // bit_length of the IMT buffer
 *  - imt_bit_length: imt // imt buffer
 * @param fmd
 * @param buf
 * @param buf_size
 */
void fmd_fmd_serialize_to_buffer(fmd_fmd_t *fmd, unsigned char **buf_ret, uint64_t *buf_size_re);
void fmd_fmd_serialize_to_buffer_imt_helper(fmd_imt_node_t *node, fmd_bs_t *bs, uint_t *widx);


#endif //FMD_FMD_FMD_H
