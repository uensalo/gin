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
 * It uses a cache-optimised bitvector:
 *      Each superblock has two blocks.
 *      Each block contains counts for 256 bits.
 *      Hence, each superblock contains 64 bytes and bitvector accesses are cache aligned.
 *      Each superblock header is 64 bytes:
 *          7 chars * 8 byte counts each + 7 chars * 1 byte count each + 1 byte carry bits
 *
 * The L column is stored as flat wavelets in six bitvectors. First bitvector is set
 *
 * Buffer write order:
 * no_chars (8 bytes)
 * inverse suffix array sampling rate (8 bytes)
 * f column cumsum (64 bytes, 8 bytes per character)
 * headers (64 bytes per superblock, N/64 superblocks)
 *      header structure 7 chars x 8 byte counts + 7 chars x 1 byte count + 1 byte "full"
 * primary bitvectors (4, non-interleaved, in ceil(N/8) bytes)
 * secondary bitvectors (2, non-interleaved, in ceil(N/8) bytes)
 * sa occupancy header (1, in 8*ceil(N/8) bytes)
 * sa occupancy bitvector (1, in ceil(N/8) bytes)
 * sa values (8 bytes per entry)
 * permutation bitvector (1, in ceil(M/8) bytes)
 */

#define DFMI_SET_BIT(ARRAY, K) ((ARRAY)[(K) >> 6] |= (UINT64_C(1) << ((K) & 63)))
#define DFMI_GET_BIT(ARRAY, K) (((ARRAY)[(K) >> 6] >> ((K) & 63)) & 1)
#define DFMI_POPCOUNT(X) __builtin_popcountll(X)

typedef union gin_dfmi_header_ {
    uint64_t buf[8];
    struct {
        uint64_t sblock_ct[7];
        uint8_t block_ct[7];
        uint8_t carry;
    } cache;
} gin_dfmi_header_t;

typedef struct gin_dfmi_ {
    // stored fields
    uint64_t size;      // size in bytes of the dfmi
    uint64_t *buf;      // buffer containing *everything*
    // non-stored fields: these are in-memory pointers for convenience, derived from the buffer
    uint64_t no_chars;  // includes the terminator and the permutation
    uint64_t isa_r;     // suffix array sampling rate
    uint64_t *f;        // 8 element f column cumsum:   $, c0, c1, A, C, G, N, T
                        // 0,1 are skipped over - still reflected in A
    gin_dfmi_header_t *headers;  // superblocks & blocks
    uint64_t *bv_1[4];  // primary bitmaps:   00x,01x,10x,11x
    uint64_t *bv_2[2];  // secondary bitmaps: xx0,xx1
    uint64_t *bv_sa_h;  // sa headers: popcount cache
    uint64_t *bv_sa_occ;// sa bitmap:  bit set if isa sampled
    uint64_t *sa;       // sa entries themselves
    uint64_t *bv_p;     // perm bitmap: stores permutation entries, needed for decoding
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
