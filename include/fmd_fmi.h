#ifndef FMD_FMD_FMI_H
#define FMD_FMD_FMI_H
#include "fmd_common.h"
#include "fmd_bitstream.h"
#include "fmd_string.h"
#include "fmd_table.h"
#include "fmd_vector.h"
#include "fmd_tree.h"
#include "divsufsort64.h"

typedef uint64_t count_t;

#define FMD_FMI_ALPHABET_SIZE_BIT_LENGTH 40
#define FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH 40
#define FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH 40
#define FMD_FMI_RNK_SAMPLE_RATE_BIT_LENGTH 40
#define FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH 40
#define FMD_FMI_CHAR_COUNT_BIT_LENGTH 40

typedef struct fmd_fmi_ {
    int_t no_chars;            // number of total characters encoded by the FMI
    int_t no_chars_per_block;  // number of characters per block
    int_t isa_sample_rate;     // sampling rate of the suffix array
    int_t no_bits_per_char;    // number of bits used per char
    int_t alphabet_size;       // the encoding will always be ceil(log2 |alphabet_size|) per char
    fmd_vector_t *alphabet;    // vector storing the characters in the alphabet in sorted order
    fmd_table_t *c2e;          // char to encoding
    fmd_table_t *e2c;          // encoding to char
    count_t *char_counts;      // the array storing cumulative character count for FMI querying
    fmd_table_t *isa;          // inverse suffix array sample
    int_t bv_start_offset;     // start of [counts]-[block_payload] chain in the bitvector
    fmd_bs_t *bits;            // buffer storing everything, including the members of this struct for serializability
} fmd_fmi_t;

typedef struct {
    int_t lo;
    int_t hi;
    int_t pos;
    fmd_string_t *pattern;
} fmd_fmi_qr_t; // query record

void fmd_fmi_init(fmd_fmi_t **fmi,
                  fmd_string_t *string,
                  int_t rank_sample_rate,
                  int_t isa_sample_rate);
count_t fmd_fmi_query_count(fmd_fmi_t *fmi, fmd_string_t *string);
fmd_vector_t *fmd_fmi_query_locate(fmd_fmi_t *fmi, fmd_string_t *string);

// helpers
typedef struct fmd_fmi_init_charset_flatten_p_ {
    fmd_vector_t *chars;
    fmd_table_t *e2c;
    fmd_table_t *c2e;
} fmd_fmi_init_charset_flatten_p_t;
void fmd_fmi_init_charset_flatten(void *key, void *value, void *params); //(*ftrav_kv)(void *key, void *value, void *p);
typedef struct fmd_fmi_init_isa_write_p_ {
    fmd_bs_t *bits;
    word_t isa_sampling_rate;
    uint_t base_bit_idx;
} fmd_fmi_init_isa_write_p_t;
void fmd_fmi_init_isa_write(void *key, void *value, void *params); //(*ftrav_kv)(void *key, void *value, void *p);
bool fmd_fmi_advance_query(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr);
count_t fmd_fmi_rank(fmd_fmi_t *fmi, word_t enc, int_t pos);
fmd_vector_t *fmd_fmi_sa(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr);
word_t fmd_fmi_get(fmd_fmi_t *fmi, int_t pos);

fmd_fmi_t* fmd_fmi_copy(fmd_fmi_t* fmi);
void fmd_fmi_free(fmd_fmi_t *fmi);
uint_t fmd_fmi_hash(fmd_fmi_t *fmi);
int fmd_fmi_comp(fmd_fmi_t *f1, fmd_fmi_t *f2);

static fmd_fstruct_t fmd_fstruct_fmi = {
    fmd_fmi_comp,
    fmd_fmi_hash,
    fmd_fmi_free,
    fmd_fmi_copy
};

#endif //FMD_FMD_FMI_H
