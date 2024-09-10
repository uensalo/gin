#ifndef GIN_DNA_FMI_H
#define GIN_DNA_FMI_H

#include <stdint.h>
#include <stdlib.h>


/*
 * gin_dna_fmi implements a cache-optimised FM index with a storage requirement of < 2N bytes.
 *
 * In particular, gin_dna_fmi is optimised for an alphabet size of 8, which is the bare minimum
 * required for gin_gin over DNA ($, c0, c1, A, C, G, T, N). The permutation is implicitly stored.
 *
 */

#define DFMI_SET_BIT(ARRAY, K) ((ARRAY)[(K) >> 6] |= (UINT64_C(1) << ((K) & 63)))
#define DFMI_GET_BIT(ARRAY, K) (((ARRAY)[(K) >> 6] >> ((K) & 63)) & 1)
#define DFMI_POPCOUNT(X) __builtin_popcountll(X)

#define DFMI_NO_CHARS 6

#define DFMI_A 0
#define DFMI_C 1
#define DFMI_G 2
#define DFMI_T 3
#define DFMI_N 4
#define DFMI_X 5

const static char gin_dfmi_dec[DFMI_NO_CHARS] = {
        'A', 'C', 'G', 'T', 'N', '('
};

const static int8_t gin_dfmi_enc[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1,  DFMI_X, -1, -1, -1, -1, -1, -1, -1,
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

inline static uint64_t gin_dfmi_wavelet_enc(uint8_t enc, const uint64_t *bv) {
    switch (enc) {
        case DFMI_A:  // A: 001
            return ((!bv[0]) & (!bv[1]) &   bv[2]);
        case DFMI_C:  // C: 010
            return ((!bv[0]) &   bv[1]  & (!bv[2]));
        case DFMI_G:  // G: 011
            return ((!bv[0]) &   bv[1]  &   bv[2]);
        case DFMI_T:  // T: 100
            return (bv[0]    & (!bv[1]) & (!bv[2]));
        case DFMI_N:  // N: 101
            return (bv[0]    & (!bv[1]) &   bv[2]);
        case DFMI_X:  // X: 110
            return ((!bv[0]) & (!bv[1]) &   bv[2]);
        default:
            fprintf(stderr, "[gin_dna_fmi]: invalid encoding\n");
            return -1;
    }
}

/*
 * 40 bit unsigned integer for popcounts
 */
typedef struct gin_uint40_ {
    uint8_t bytes[5];
} gin_uint40_t;

inline static uint64_t gin_uint40_get(gin_uint40_t u) {
    uint64_t result = 0;
    memcpy(&result, u.bytes, sizeof(u.bytes));
    return result;
}

inline static void gin_uint40_set(gin_uint40_t *u, uint64_t val) {
    memcpy(u->bytes, &val, sizeof(u->bytes));
}

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
 *  - X: 110: !bv[0] &  !bv[1] &  bv[2]
 */
typedef union gin_dfmi_sbc_ {
    uint64_t buf[16]; // 128 bytes per superblock
    struct {
        gin_uint40_t sbc[DFMI_NO_CHARS]; // 6x5 byte integers, can index up to 1 TB
        uint8_t pad[2];                // NOT used, here for cache alignment
        struct {
            uint8_t  bc[DFMI_NO_CHARS];  // cumulative sum of each char in the block
            uint8_t  pad[2]; // NOT used, here for cache alignment
            uint64_t bv[3];  // bitvectors for rank access
        } blocks[3];
    } sbc;
} gin_dfmi_sbc_t;

/*
 * Cache aligned suffix array occupancy bitvector - fits in one cache line.
 * popcnt stores the number of set bits up until the beginning of the block.
 */
typedef union gin_dfmi_sao_ {
    uint64_t buf[8];
    struct {
        gin_uint40_t popcnt; // popcount until this point
        uint16_t bc[5];      // cumulative in-block popcount
        uint64_t bv[6];      // actual bitvectors
    } sa_occ_sbc;
} gin_dfmi_sao_t;

/*
 * The L column is stored as flat wavelets in three bitvectors.
 *
 * Buffer write order:
 * no_chars (8 bytes)
 * inverse suffix array sampling rate (8 bytes)
 * Position of original string terminator in the BWT (8 bytes)
 * 40 byte padding for cache alignment, may be used for future use
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
 * permutation bitvector (1, in ceil(M/8) bytes)
 */
typedef struct gin_dfmi_ {
    // stored fields
    uint64_t size;          // size in bytes of the dfmi
    uint64_t *buf;          // buffer containing *everything*
    // non-stored fields: these are in-memory pointers for convenience, derived from the buffer
    uint64_t no_chars;      // includes the terminator and the permutation
    uint64_t isa_r;         // suffix array sampling rate
    uint64_t str_term_pos;  // string terminator position in L
    uint64_t *f;            // 8 element f column cumsum:   $, c0, c1, A, C, G, T, N
    gin_dfmi_sbc_t *lf;     // superblocks & blocks
    gin_dfmi_sao_t *sa_occ; // sa headers: popcount cache
    uint64_t *sa;           // sa entries themselves
    uint64_t *bv_p;         // perm bitmap: stores permutation entries, needed for decoding
} gin_dfmi_t;

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
