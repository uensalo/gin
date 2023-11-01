/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_fmi.h is part of gin
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
#ifndef GIN_GIN_FMI_H
#define GIN_GIN_FMI_H
#include "gin_common.h"
#include "gin_bitstream.h"
#include "gin_string.h"
#include "gin_table.h"
#include "gin_vector.h"
#include "gin_tree.h"
#include "divsufsort64.h"

typedef uint64_t count_t;

#define GIN_FMI_MAX_ALPHABET_SIZE 256
#define GIN_FMI_NO_BITS_BIT_LENGTH 64
#define GIN_FMI_ALPHABET_SIZE_BIT_LENGTH 40
#define GIN_FMI_ALPHABET_ENTRY_BIT_LENGTH 40
#define GIN_FMI_ALPHABET_ENCODING_BIT_LENGTH 40
#define GIN_FMI_RNK_SAMPLE_RATE_BIT_LENGTH 40
#define GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH 64 // 64 is better than 40 due to word alignment
#define GIN_FMI_CHAR_COUNT_BIT_LENGTH 40

#define GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH 64 // has to be word aligned
#define GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH 448 // has to be word aligned
#define GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE 9

typedef struct gin_fmi_ {
    int_t no_chars;            // number of total characters encoded by the FMI
    int_t no_chars_per_block;  // number of characters per block
    int_t isa_sample_rate;     // sampling rate of the suffix array
    int_t no_bits_per_char;    // number of bits used per char
    int_t alphabet_size;       // the encoding will always be ceil(log2 |alphabet_size|) per char
    int_t *alphabet;           // vector storing the characters in the alphabet in sorted order
    int_t *c2e;                // char to encoding
    int_t *e2c;                // encoding to char
    count_t *char_counts;      // the array storing cumulative character count for FMI querying
    gin_bs_t *bits;            // buffer storing everything, including the members of this struct for serializability
    // not stored in bits, helper fields
    int_t sa_start_offset;     // number of bits into the bitstream indicating the starting index of the suffix array entries
    int_t sa_bv_start_offset;  // number of bits into the bitstream indicating the starting index of the SA occupancy bv
    int_t bv_start_offset;     // start of [counts]-[block_payload] chain in the bitvector
    int_t no_bits;             // number of bits written to the bitstream, possibly used for serialization
} gin_fmi_t;

typedef struct {
    int_t lo;
    int_t hi;
    int_t pos;
    gin_string_t *pattern;
} gin_fmi_qr_t; // query record
gin_fmi_qr_t *gin_fmi_qr_init(int_t lo, int_t hi, int_t pos, gin_string_t *pattern);
void gin_fmi_qr_free(gin_fmi_qr_t *qr);

void gin_fmi_init(gin_fmi_t **fmi,
                  gin_string_t *string,
                  int_t rank_sample_rate,
                  int_t isa_sample_rate);

void gin_fmi_init_with_sa(gin_fmi_t **fmi,
                          gin_string_t *string,
                          int_t *sa,
                          int_t rank_sample_rate,
                          int_t isa_sample_rate);
count_t gin_fmi_query_count(gin_fmi_t *fmi, gin_string_t *string);
gin_vector_t *gin_fmi_query_locate(gin_fmi_t *fmi, gin_string_t *string);

typedef struct gin_fmi_init_isa_write_p_ {
    gin_bs_t *bits;
    word_t isa_sampling_rate;
    int_t sa_start_offset;
    int_t sa_bv_start_offset;
    int_t no_traversed;
} gin_fmi_init_isa_write_p_t;
void gin_fmi_init_isa_write(void *key, void *value, void *params); //(*ftrav_kv)(void *key, void *value, void *p);
bool gin_fmi_advance_query(gin_fmi_t *fmi, gin_fmi_qr_t *qr);
bool gin_fmi_query_precedence_range(gin_fmi_t *fmi, gin_fmi_qr_t *qr, char_t c, int_t *lo, int_t *hi);
count_t gin_fmi_rank(gin_fmi_t *fmi, word_t enc, int_t pos);
gin_vector_t *gin_fmi_sa(gin_fmi_t *fmi, gin_fmi_qr_t *qr);
word_t gin_fmi_get(gin_fmi_t *fmi, int_t pos);

gin_fmi_t* gin_fmi_copy(gin_fmi_t* fmi);
void gin_fmi_free(gin_fmi_t *fmi);
void gin_fmi_free_disown(gin_fmi_t *fmi);
uint_t gin_fmi_hash(gin_fmi_t *fmi);
int gin_fmi_comp(gin_fmi_t *f1, gin_fmi_t *f2);

static gin_fstruct_t gin_fstruct_fmi = {
        (fcomp)gin_fmi_comp,
        (fhash)gin_fmi_hash,
        (ffree)gin_fmi_free,
        (fcopy)gin_fmi_copy
};

// sdsl-like interface
void gin_fmi_serialize_from_buffer(unsigned char *buf, uint64_t buf_size, gin_fmi_t **gin_ret);
void gin_fmi_serialize_to_buffer(gin_fmi_t *gin, unsigned char **buf_ret, uint64_t *buf_size_re);
void gin_fmi_bwt(gin_fmi_t *fmi, uint64_t *buf, uint64_t start, uint64_t end);

#endif //GIN_GIN_FMI_H
