#ifndef GIN_DNA_FMI_H
#define GIN_DNA_FMI_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


/*
 * gin_dna_fmi implements a cache-optimised FM index with a storage requirement of < 2N bytes.
 *
 * In particular, gin_dna_fmi is optimised for an alphabet size of 8, which is the bare minimum
 * required for gin_gin over DNA ($, c0, c1, A, C, G, T, N). The permutation is implicitly stored.
 *
 */

#define DFMI_POPCOUNT(X) __builtin_popcountll(X)
#define DFMI_POPCOUNT_MASK(m) ((UINT64_C(1) << ((m) + 1)) - 1)

#define DFMI_NO_RANKED_CHARS 6
#define DFMI_NO_ALL_CHARS 10

#define DFMI_X 0
#define DFMI_A 1
#define DFMI_C 2
#define DFMI_G 3
#define DFMI_N 4
#define DFMI_T 5

// Non-ranked characters
#define DFMI_Y 6
#define DFMI_$ 7
#define DFMI_a 8
#define DFMI_b 9

#define DFMI_CEIL(a,b) (((a)-1) / (b) + 1)

const static char gin_dfmi_dec_lex[DFMI_NO_ALL_CHARS] = {
        '\0', '(', ')', ',', '.', 'A', 'C', 'G', 'N', 'T'
};

const static char gin_dfmi_dec[DFMI_NO_ALL_CHARS] = {
        '(', 'A', 'C', 'G', 'N', 'T', ')', '\0', ',', '.'
};

const static uint8_t gin_dfmi_enc[256] = {
        DFMI_$, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1,  DFMI_X, DFMI_Y, -1, -1, DFMI_a, -1, DFMI_b, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1,  DFMI_A, -1,  DFMI_C, -1, -1, -1,  DFMI_G, -1, -1, -1, -1, -1, -1,  DFMI_N, -1,
        -1, -1, -1, -1,  DFMI_T, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

/*
 * 40 bit unsigned integer for popcounts
 */
typedef struct gin_uint40_ {
    uint8_t bytes[5];
} gin_uint40_t;

#define DFMI_UINT40_GET(U) ({ \
    uint64_t result = 0; \
    memcpy(&result, (U).bytes, sizeof((U).bytes)); \
    result; \
})

#define DFMI_UINT40_SET(U, VAL) \
    memcpy((U)->bytes, &(VAL), sizeof((U)->bytes))

typedef union gin_dfmi_header_ {
    uint64_t buf[8];
    struct {
        uint64_t fmi_size_in_bytes;
        uint64_t no_chars;
        uint64_t isa_rate;
        uint64_t str_term_pos;
        uint64_t no_dfmi_x;
        uint64_t no_sa_values;
        uint64_t padding[2]; // reserved
    } values;
} gin_dfmi_header_t;

/*
 * Cache aligned superblock - fits in two cache lines (2x64).
 * Can prefetch the second line before a memory access is made if necessary
 * Stores counts of 6 characters: A, C, G, T, N, c0. Other chars are not used
 * for rank queries. The ranks are computed as follows: Each char is assigned
 * a codeword, which indicates how the char is to be recovered from the bvs:
 *  - A: 001: !bv[0] &  !bv[1] &  bv[2]
 *  - C: 010: !bv[0] &   bv[1] & !bv[2]
 *  - G: 011: !bv[0] &   bv[1] &  bv[2]
 *  - T: 100:  bv[0] &  !bv[1] & !bv[2]
 *  - N: 101:  bv[0] &  !bv[1] &  bv[2]
 *  - X: 110: bv[0] &    bv[1] &  bv[2]
 */
typedef union gin_dfmi_sb_ {
    uint64_t buf[16]; // 128 bytes per superblock
    struct {
        gin_uint40_t sbc[DFMI_NO_RANKED_CHARS]; // 6x5 byte integers, can index up to 1 TB
        uint8_t pad[2];                         // NOT used, here for cache alignment
        struct {
            uint8_t  bc[DFMI_NO_RANKED_CHARS];  // cumulative sum of each char in the block
            uint8_t  pad[2]; // NOT used, here for cache alignment
            uint64_t bv[3];  // bitvectors for rank access
        } blocks[3];
    } sb;
} gin_dfmi_sb_t;

inline static uint64_t gin_dfmi_wavelet_enc(uint8_t enc, const uint64_t *bv) {
    switch (enc) {
        case DFMI_X:  // X: 001
            return ((!bv[0]) & (!bv[1]) & bv[2]);
        case DFMI_A:  // A: 010
            return ((!bv[0]) & bv[1] & (!bv[2]));
        case DFMI_C:  // C: 011
            return ((!bv[0]) & bv[1] & bv[2]);
        case DFMI_G:  // G: 100
            return (bv[0] & (!bv[1]) & (!bv[2]));
        case DFMI_N:  // N: 101
            return (bv[0] & (!bv[1]) & bv[2]);
        case DFMI_T:  // T: 110
            return (bv[0] & bv[1] & (!bv[2]));
        default:
            fprintf(stderr, "[gin_dna_fmi.h]: invalid encoding\n");
            return -1;
    }
}

#define DFMI_L_SET_BIT(SB, CH, K) do { \
    int sb_index = (K) / 192; \
    int block_index = ((K) % 192) / 64; \
    int block_offset = ((K) % 192) % 64; \
    switch (CH) { \
        case DFMI_A: \
            (SB)[sb_index].sb.blocks[block_index].bv[2] |= (UINT64_C(1) << block_offset); \
            break; \
        case DFMI_C: \
            (SB)[sb_index].sb.blocks[block_index].bv[1] |= (UINT64_C(1) << block_offset); \
            break; \
        case DFMI_G: \
            (SB)[sb_index].sb.blocks[block_index].bv[1] |= (UINT64_C(1) << block_offset); \
            (SB)[sb_index].sb.blocks[block_index].bv[2] |= (UINT64_C(1) << block_offset); \
            break; \
        case DFMI_T: \
            (SB)[sb_index].sb.blocks[block_index].bv[0] |= (UINT64_C(1) << block_offset); \
            break; \
        case DFMI_N: \
            (SB)[sb_index].sb.blocks[block_index].bv[0] |= (UINT64_C(1) << block_offset); \
            (SB)[sb_index].sb.blocks[block_index].bv[2] |= (UINT64_C(1) << block_offset); \
            break; \
        case DFMI_X: \
            (SB)[sb_index].sb.blocks[block_index].bv[0] |= (UINT64_C(1) << block_offset); \
            (SB)[sb_index].sb.blocks[block_index].bv[1] |= (UINT64_C(1) << block_offset); \
            break; \
    } \
} while (0)

#include <stdint.h>

#define DFMI_L_CHAR(SB, j) ( \
    ((((SB)[(j) / 192].sb.blocks[((j) % 192) / 64].bv[0] >> ((j) % 64)) & UINT64_C(1)) << 2) + \
    ((((SB)[(j) / 192].sb.blocks[((j) % 192) / 64].bv[1] >> ((j) % 64)) & UINT64_C(1)) << 1) + \
    ((((SB)[(j) / 192].sb.blocks[((j) % 192) / 64].bv[2] >> ((j) % 64)) & UINT64_C(1))) - 1 \
)


/*
 * Cache aligned suffix array occupancy bitvector - fits in one cache line.
 * popcnt stores the number of set bits up until the beginning of the block.
 */
typedef union gin_dfmi_sao_ {
    uint64_t buf[8];
    struct {
        gin_uint40_t popcnt;  // popcount until this point
        uint16_t bc[5];       // cumulative in-block popcount
        uint64_t bv[6];       // actual bitvectors
        uint8_t pad;
    } sa_occ_sbc;
} gin_dfmi_sao_t;

#define DFMI_SAO_SET_BIT(SAO, K) \
    ((SAO)[(K) / 384].sa_occ_sbc.bv[((K) % 384) / 64] |= (UINT64_C(1) << (((K) % 384) % 64)))

#define DFMI_SAO_GET_BIT(SAO, K) \
    ((((SAO)[(K) / 384].sa_occ_sbc.bv[((K) % 384) / 64] >> (((K) % 384) % 64)) & UINT64_C(1)))

#define DFMI_SAO_RANK(SAO, j) ( \
    { \
        uint64_t s = (j) / 384; \
        uint64_t b = ((j) % 384) / 64; \
        uint64_t m = (j) % 64; \
        uint64_t popcnt = DFMI_UINT40_GET((SAO)[s].sa_occ_sbc.popcnt); \
        uint64_t block_cache = b ? (uint64_t)(SAO)[s].sa_occ_sbc.bc[b - 1] : 0; \
        uint64_t bv_popcount = DFMI_POPCOUNT((SAO)[s].sa_occ_sbc.bv[b] & DFMI_POPCOUNT_MASK(m)); \
        popcnt + block_cache + bv_popcount; \
    } \
)

/*
 * The L column is stored as flat wavelets in three bitvectors.
 *
 * Buffer write order:
 * FMI Header (64 bytes, see gin_dfmi_header)
 * F column cumsum (64 bytes, 8 bytes per character)
 * L column (2 cache lines per superblock):
 *      128-byte superblocks, ceil(N/192) many:
 *          32 byte sblock header:
 *              6x5 byte cumulative counts
 *              2 byte padding
 *          3 x 32 byte blocks:
 *              6 byte cumulative counts
 *              2 byte padding
 *              24 byte bitvectors (bv0 bv1 bv2)
 * SA occupancy (1 cache line per superblock):
 *      64-byte superblocks, ceil(N/384) many:
 *          5 byte popcount
 *          4x2bytes intermediate cumulative sum
 *          6 blocks:
 *              8 byte payload
 * SA values (5 bytes per entry, ceil(N/isa_r) many)
 * // not implemented yet permutation bitvector (1, in ceil(M/8) bytes)
 */
typedef struct gin_dfmi_ {
    // stored fields
    uint64_t *buf;          // buffer containing *everything*
    // non-stored fields: these are in-memory pointers for convenience, derived from the buffer
    gin_dfmi_header_t header;
    uint64_t *f;            // 8 element f column cumsum:   $, c0, c1, A, C, G, T, N
    gin_dfmi_sb_t *l;       // superblocks & blocks
    gin_dfmi_sao_t *sa_occ; // sa headers: popcount cache
    gin_uint40_t *sa;       // sa entries, sampled
    uint64_t *unr;          // unranked bitmap: stores permutation chars, terminator, and vertex start markers
} gin_dfmi_t;

// Generic functions
gin_dfmi_t* gin_dfmi_build(const char* str, uint64_t size, uint64_t isa_rate);
int64_t gin_dfmi_rank(gin_dfmi_t *dfmi, uint64_t pos, char c);
void gin_dfmi_sa(gin_dfmi_t *dfmi, uint64_t *buf, uint64_t start, uint64_t end);
void gin_dfmi_bwt(gin_dfmi_t *dfmi, uint64_t *buf, uint64_t start, uint64_t end);
void gin_dfmi_to_buffer(gin_dfmi_t *dfmi, uint8_t **data, uint64_t *size);
void* gin_dfmi_from_buffer(uint8_t *data, uint64_t size);
int64_t gin_dfmi_bwt_length(gin_dfmi_t *dfmi);
void gin_dfmi_populate_alphabet(gin_dfmi_t *dfmi, int64_t** alphabet, int64_t *alphabet_size);
int64_t gin_dfmi_char_sa_base(gin_dfmi_t *dfmi, char c);
int64_t gin_dfmi_size_in_bytes(gin_dfmi_t *dfmi);
void gin_dfmi_decode(gin_dfmi_t *dfmi, char **string, uint64_t *len);

int gin_dfmi_comp(gin_dfmi_t *dfmi1, gin_dfmi_t *dfmi2);
uint64_t gin_dfmi_hash(gin_dfmi_t *dfmi);
void gin_dfmi_free(gin_dfmi_t *dfmi);
void* gin_dfmi_copy(gin_dfmi_t *dfmi);

#endif //GIN_DNA_FMI_H
