/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * gin_fmi.c is part of gin
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
#include "gin_fmi.h"
#include "string.h"

gin_fmi_qr_t *gin_fmi_qr_init(int_t lo, int_t hi, int_t pos, gin_string_t *pattern) {
    gin_fmi_qr_t *qr = calloc(1, sizeof(gin_fmi_qr_t));
    qr->lo = lo;
    qr->hi = hi;
    qr->pos = pos;
    qr->pattern = pattern;
    return qr;
}
void gin_fmi_qr_free(gin_fmi_qr_t *qr) {
    if(qr) free(qr);
}

void gin_fmi_init(gin_fmi_t **fmi,
                  gin_string_t *string,
                  int_t rank_sample_rate,
                  int_t isa_sample_rate) {
    int64_t *sa = calloc(string->size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)string->seq, (saidx64_t*)sa, string->size+1);
    gin_fmi_init_with_sa(fmi,string,sa,rank_sample_rate,isa_sample_rate);
    free(sa);
}

void gin_fmi_init_with_sa(gin_fmi_t **fmi,
                          gin_string_t *string,
                          int_t *sa,
                          int_t rank_sample_rate,
                          int_t isa_sample_rate) {
    gin_fmi_t *f = calloc(1, sizeof(gin_fmi_t));
    if(!f) {
        *fmi = NULL;
        return;
    }
    f->no_chars_per_block = rank_sample_rate;
    f->isa_sample_rate = isa_sample_rate;
    f->no_chars = string->size + 1; // +1 due to terminator
    /**************************************************************************
     * Step 1 - Construct a dictionary for characters and character counts
    **************************************************************************/
    /******************************************************
    * Step 1a - Tabulate unique characters
    ******************************************************/
    char occ[GIN_FMI_MAX_ALPHABET_SIZE];
    memset(occ, 0, GIN_FMI_MAX_ALPHABET_SIZE);
    for(int_t i = 0; i < f->no_chars; i++) {
        uint64_t c = (uint64_t)string->seq[i];
        occ[c] = 1;
    }
    /******************************************************
    * Step 1b - Fill encoding and decoding tables
    ******************************************************/
    int_t alphabet_size = 0;
    for(int_t i = 0; i < GIN_FMI_MAX_ALPHABET_SIZE; i++) {
        alphabet_size += occ[i];
    }
    f->alphabet_size = alphabet_size;
    f->no_bits_per_char = gin_ceil_log2(f->alphabet_size);
    f->alphabet = calloc(alphabet_size, sizeof(int_t));
    f->e2c = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    f->c2e = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    int_t t = 0;
    for(int_t i = 0; i < GIN_FMI_MAX_ALPHABET_SIZE; i++) {
        if(occ[i]) {
            f->alphabet[t] = i;
            f->c2e[i] = t;
            f->e2c[t] = i;
            t++;
        }
    }
    /**************************************************************************
    * Step 2 - Compute the suffix array and the BWT of the input and sample
    *************************************************************************/
    gin_tree_t *sa_tree;
    gin_tree_init(&sa_tree, &prm_fstruct, &prm_fstruct); // tree sort
    if(!sa_tree) {
        free(f->alphabet);
        free(f->c2e);
        free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    gin_string_t *bwt;
    gin_string_init(&bwt, (int_t)f->no_chars);
    for(int_t i = 0; i < (int_t)f->no_chars; i++) {
        int64_t sa_val = sa[i];
        bwt->seq[i] = sa_val ? string->seq[sa_val-1] : GIN_STRING_TERMINATOR;
        if((sa_val % f->isa_sample_rate) == 0) {
            gin_tree_insert(sa_tree, (void*)i, (void*)sa_val);
        }
    }
    bwt->size = string->size+1;
    /**************************************************************************
    * Step 3 - Encode everything into the bitvector (blocks, ranks, sa, abc)
    *************************************************************************/
    gin_bs_t *bits;
    gin_bs_init(&bits);
    if(!bits) {
        gin_string_free(bwt);
        free(f->alphabet);
        free(f->c2e);
        free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    /******************************************************
    * Step 3a - Write the header
    ******************************************************/
    int_t widx = GIN_FMI_NO_BITS_BIT_LENGTH; // skip the bits field
    gin_bs_write_word(bits, widx, f->no_chars, GIN_FMI_CHAR_COUNT_BIT_LENGTH);
    widx+=GIN_FMI_CHAR_COUNT_BIT_LENGTH;
    gin_bs_write_word(bits, widx, f->no_chars_per_block, GIN_FMI_RNK_SAMPLE_RATE_BIT_LENGTH);
    widx+=GIN_FMI_RNK_SAMPLE_RATE_BIT_LENGTH;
    gin_bs_write_word(bits, widx, f->isa_sample_rate, GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH);
    widx+=GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;
    gin_bs_write_word(bits, widx, f->alphabet_size, GIN_FMI_ALPHABET_SIZE_BIT_LENGTH);
    widx+=GIN_FMI_ALPHABET_SIZE_BIT_LENGTH;
    /******************************************************
    * Step 3b - Write the alphabet
    ******************************************************/
    for(int_t i = 0; i < f->alphabet_size; i++) {
        gin_bs_write_word(bits, widx, (word_t)f->alphabet[i], GIN_FMI_ALPHABET_ENTRY_BIT_LENGTH);
        widx+=GIN_FMI_ALPHABET_ENTRY_BIT_LENGTH;
        word_t encoding = f->c2e[f->alphabet[i]];
        gin_bs_write_word(bits, widx, encoding, GIN_FMI_ALPHABET_ENCODING_BIT_LENGTH);
        widx+=GIN_FMI_ALPHABET_ENCODING_BIT_LENGTH;
    }
    /****************************
    * Align to word boundary
    ****************************/
    widx= (1 + ((widx - 1) >> WORD_LOG_BITS)) << WORD_LOG_BITS;
    /******************************************************
    * Step 3c - Write the sampled suffix array entries
    ******************************************************/
    /******************************************************
    * Step 3d - Write the suffix array occupancy bitvector
    ******************************************************/
    // first write the suffix array entries
    f->sa_start_offset = widx;
    widx += GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH*(1+f->no_chars/f->isa_sample_rate);
    f->sa_bv_start_offset = widx;
    int_t sa_bv_no_blocks = 1 + (f->no_chars - 1) / GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH;
    // reserve space for the bitvector with a dud write
    gin_bs_write_word(bits, f->sa_bv_start_offset + (sa_bv_no_blocks << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE), 0, 1);
    gin_fmi_init_isa_write_p_t write_p;
    write_p.bits = bits;
    write_p.sa_start_offset = f->sa_start_offset;
    write_p.sa_bv_start_offset = f->sa_bv_start_offset;
    write_p.isa_sampling_rate = f->isa_sample_rate;
    write_p.no_traversed = 0;
    gin_tree_inorder(sa_tree, &write_p, gin_fmi_init_isa_write);
    gin_tree_free(sa_tree);

    // write the popcounts
    int_t sa_ridx = f->sa_bv_start_offset;
    uint_t sa_bv_cur_popcnt = 0;
    int_t no_words_per_block = GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH >> WORD_LOG_BITS;
    for(int_t i = 0; i < sa_bv_no_blocks; i++) {
        gin_bs_write_word(bits, sa_ridx, sa_bv_cur_popcnt, GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH);
        sa_ridx += GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH;
        for(int_t j = 0; j < no_words_per_block; j++) {
            word_t payload = 0;
            gin_bs_read_word(bits, sa_ridx, WORD_NUM_BITS, &payload);
            sa_ridx += WORD_NUM_BITS;
            sa_bv_cur_popcnt += gin_popcount64(payload);
        }
    }
    int_t sa_bv_occ_size = sa_bv_no_blocks << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE;
    widx+=sa_bv_occ_size;
    /****************************
    * Align to word boundary
    ****************************/
    widx=(1+((widx-1)>>WORD_LOG_BITS))<<WORD_LOG_BITS;
    /******************************************************
    * Step 3e - Write rank caches and payloads
    ******************************************************/
    f->bv_start_offset = widx;
    int_t no_blocks = 1 + (bwt->size - 1) / f->no_chars_per_block;
    count_t *block_count = calloc(f->alphabet_size, sizeof(count_t));
    for(int_t c = 0; c < f->alphabet_size; c++) {
        gin_bs_write_word(bits, widx, block_count[c], GIN_FMI_CHAR_COUNT_BIT_LENGTH);
        widx += GIN_FMI_CHAR_COUNT_BIT_LENGTH;
    }
    for(int_t i = 0; i < no_blocks; i++) {
        int_t s = i * (int_t)f->no_chars_per_block;
        int_t e = MIN2((i + 1) * (int_t)f->no_chars_per_block,bwt->size);
        for(int_t j = s; j < e; j++) {
            int_t enc = f->c2e[bwt->seq[j]];
            ++block_count[enc];
            gin_bs_write_word(bits, widx, enc, f->no_bits_per_char);
            widx += f->no_bits_per_char;
        }
        for(int_t c = 0; c < f->alphabet_size; c++) {
            gin_bs_write_word(bits, widx, block_count[c], GIN_FMI_CHAR_COUNT_BIT_LENGTH);
            widx += GIN_FMI_CHAR_COUNT_BIT_LENGTH;
        }
    }
    /******************************************************
    * Step 3f - Compute a cumulative sum of char counts
    ******************************************************/
    count_t cum = 0;
    f->char_counts = calloc(f->alphabet_size, sizeof(count_t));
    for(int_t i = 0; i < f->alphabet_size; i++) {
        f->char_counts[i] = cum;
        cum += block_count[i];
    }
    free(block_count);
    gin_bs_fit(bits, widx);
    f->bits = bits;
    f->no_bits = widx;
    gin_bs_write_word(bits, 0, widx, GIN_FMI_NO_BITS_BIT_LENGTH);
    gin_string_free(bwt);
    *fmi = f;
}

count_t gin_fmi_query_count(gin_fmi_t *fmi, gin_string_t *string) {
    gin_fmi_qr_t query;
    query.lo = 0;
    query.hi = fmi->no_chars+1;
    query.pos = string->size-1;
    query.pattern = string;
    while(query.pos > -1) {
        bool non_empty = gin_fmi_advance_query(fmi, &query);
        if(!non_empty) return 0;
    }
    return query.hi - query.lo;
}

gin_vector_t *gin_fmi_query_locate(gin_fmi_t *fmi, gin_string_t *string) {
    gin_fmi_qr_t query;
    query.lo = 0;
    query.hi = fmi->no_chars+1;
    query.pos = string->size-1;
    query.pattern = string;
    while(query.pos > -1) {
        bool non_empty = gin_fmi_advance_query(fmi, &query);
        if(!non_empty) break;
    }
    return gin_fmi_sa(fmi, &query);
}

void gin_fmi_init_isa_write(void *key, void *value, void *params) {
    // key: rank, value: offset
    int_t rank = (int_t)key;
    int_t offset = (int_t)value;
    gin_fmi_init_isa_write_p_t *p = (gin_fmi_init_isa_write_p_t*)params;
    uint_t offset_widx = p->no_traversed * GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH + p->sa_start_offset;
    gin_bs_write_word(p->bits, offset_widx, offset, GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH);

    int_t div = rank / GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH; // div is the index of the block
    int_t rem = rank % GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH;
    uint_t sa_occ_widx = p->sa_bv_start_offset + GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH + (div << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE) + rem;
    gin_bs_write_word(p->bits, sa_occ_widx, (word_t)1, 1);
    p->no_traversed++;
}

bool gin_fmi_advance_query(gin_fmi_t *fmi, gin_fmi_qr_t *qr) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    encoding = fmi->c2e[qr->pattern->seq[qr->pos]];
    count_t rank_lo_m_1 = qr->lo ? gin_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? gin_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    qr->lo = (int_t)(base + rank_lo_m_1);
    qr->hi = (int_t)(base + rank_hi_m_1);
    --qr->pos;
    return qr->hi > qr->lo;
}

bool gin_fmi_query_precedence_range(gin_fmi_t *fmi, gin_fmi_qr_t *qr, char_t c, int_t *lo, int_t *hi) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    encoding = fmi->c2e[qr->pattern->seq[qr->pos]];
    count_t rank_lo_m_1 = qr->lo ? gin_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? gin_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    *lo = (int_t)(base + rank_lo_m_1);
    *hi = (int_t)(base + rank_hi_m_1);
    return *hi > *lo;
}

count_t gin_fmi_rank(gin_fmi_t *fmi, word_t enc, int_t pos) {
    // todo: read from short side
    uint_t block_idx = pos / fmi->no_chars_per_block;
    uint_t block_size = (fmi->alphabet_size * GIN_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    uint_t rank_cache_idx = fmi->bv_start_offset + block_idx * block_size;
    uint_t chars_idx = rank_cache_idx + fmi->alphabet_size * GIN_FMI_CHAR_COUNT_BIT_LENGTH;
    uint_t char_offset_into_block = pos % fmi->no_chars_per_block;
    word_t rank = 0;
    gin_bs_read_word(fmi->bits, rank_cache_idx + enc * GIN_FMI_CHAR_COUNT_BIT_LENGTH, GIN_FMI_CHAR_COUNT_BIT_LENGTH, &rank);
    for(int_t i = 0; i <= char_offset_into_block; i++) {
        word_t read_enc = 0;
        gin_bs_read_word(fmi->bits, chars_idx, fmi->no_bits_per_char, &read_enc);
        rank += enc == read_enc;
        chars_idx += fmi->no_bits_per_char;
    }
    return rank;
}

gin_vector_t *gin_fmi_sa(gin_fmi_t *fmi, gin_fmi_qr_t *qr) {
    gin_vector_t *rval;
    gin_vector_init(&rval, (int_t)qr->hi - (int_t)qr->lo, &prm_fstruct);
    if(!rval) return NULL;
    #pragma omp parallel for default(none) shared(qr, fmi, rval)
    for(int_t j = (int_t)qr->lo; j < (int_t)qr->hi; j++) {
        count_t count = 0;
        int_t i = j;
        while(1) {
            // compute physical index of the SA entry in the bv
            int_t div = i / (GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH);
            int_t rem = i % (GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH);
            int_t sa_idx = fmi->sa_bv_start_offset + GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH + (div << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE) + rem;
            word_t sampled = 0;
            gin_bs_read_word(fmi->bits, sa_idx, 1, &sampled);
            if(!sampled) {
                // sa sample not found, traverse bwt
                count++;
                // get rank of char at index i of the bwt
                word_t encoding = gin_fmi_get(fmi, (int_t)i);
                // get rank of the char at the same index
                uint64_t rank = gin_fmi_rank(fmi, encoding, (int_t)i);
                // compute the next index
                i = (int_t)(fmi->char_counts[encoding] + rank - 1);
            } else {
                // sa entry is present, interpolate and compute entry in next range
                word_t popcnt = 0;
                word_t _ = 0;
                word_t sa_entry = 0;
                // read the popcount of the sa occupancy bitvector at the position
                int_t block_base = fmi->sa_bv_start_offset + (div << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE);
                gin_bs_read_word(fmi->bits,
                                 block_base,
                                 GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH,
                                 &popcnt);
                int_t no_words_before_pos = (rem+1) >> WORD_LOG_BITS;
                int_t no_slack_bits = (rem+1)&(int_t)WORD_LOG_MASK;// ? rem & (int_t)WORD_LOG_MASK + 1 : 0;
                int_t in_block_popcnt_ridx = block_base + GIN_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH;
                for(int_t k = 0; k < no_words_before_pos; k++) {
                    _ = 0;
                    gin_bs_read_word(fmi->bits,
                                     in_block_popcnt_ridx,
                                     WORD_NUM_BITS,
                                     &_);
                    popcnt += gin_popcount64(_);
                    in_block_popcnt_ridx += WORD_NUM_BITS;
                }
                _ = 0;
                gin_bs_read_word(fmi->bits,
                                 in_block_popcnt_ridx,
                                 no_slack_bits,
                                 &_);
                popcnt += gin_popcount64(_);
                // read the SA entry at index popcnt
                gin_bs_read_word(fmi->bits,
                                 fmi->sa_start_offset + (popcnt-1)*GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH,
                                 GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH,
                                 &sa_entry);

                rval->data[j - qr->lo] = (void*)(sa_entry + count);
                break;
            }
        }
    }
    rval->size = (int_t)qr->hi - (int_t)qr->lo;
    return rval;
}

word_t gin_fmi_get(gin_fmi_t *fmi, int_t pos) {
    // todo: read from short side
    uint_t block_idx = pos / fmi->no_chars_per_block;
    uint_t block_size = (fmi->alphabet_size * GIN_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    uint_t char_idx = fmi->bv_start_offset + block_idx * block_size + fmi->alphabet_size * GIN_FMI_CHAR_COUNT_BIT_LENGTH + (pos % fmi->no_chars_per_block) * fmi->no_bits_per_char;
    word_t read_enc = 0;
    gin_bs_read_word(fmi->bits, char_idx, fmi->no_bits_per_char, &read_enc);
    return read_enc;
}

gin_fmi_t* gin_fmi_copy(gin_fmi_t* fmi) {
    gin_fmi_t *f = calloc(1, sizeof(gin_fmi_t));
    f->no_chars = fmi->no_chars;
    f->no_chars_per_block = fmi->no_chars_per_block;
    f->isa_sample_rate = fmi->isa_sample_rate;
    f->no_bits_per_char = fmi->no_bits_per_char;
    f->alphabet_size = fmi->alphabet_size;
    f->alphabet = calloc(fmi->alphabet_size, sizeof(int_t));
    memcpy(f->alphabet, fmi->alphabet, sizeof(int_t) * fmi->alphabet_size);
    f->c2e = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    memcpy(f->c2e, fmi->c2e, sizeof(int_t) * GIN_FMI_MAX_ALPHABET_SIZE);
    f->e2c = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    memcpy(f->e2c, fmi->e2c, sizeof(int_t) * GIN_FMI_MAX_ALPHABET_SIZE);
    f->char_counts = calloc(f->alphabet_size, sizeof(count_t));
    memcpy(f->char_counts, fmi->char_counts, sizeof(count_t)*f->alphabet_size);
    f->bv_start_offset = fmi->bv_start_offset;
    f->bits = gin_bs_copy(fmi->bits);
    f->no_bits = fmi->no_bits;
    return f;
}

void gin_fmi_free(gin_fmi_t *fmi) {
    if(fmi) {
        free(fmi->alphabet);
        free(fmi->c2e);
        free(fmi->e2c);
        free(fmi->char_counts);
        gin_bs_free(fmi->bits);
        free(fmi);
    }
}

void gin_fmi_free_disown(gin_fmi_t *fmi) {
    if(fmi) {
        free(fmi->alphabet);
        free(fmi->c2e);
        free(fmi->e2c);
        free(fmi->char_counts);
        gin_bs_free_disown(fmi->bits);
        free(fmi);
    }
}

uint_t gin_fmi_hash(gin_fmi_t *fmi) {
    return gin_bs_hash(fmi->bits);
}

int gin_fmi_comp(gin_fmi_t *f1, gin_fmi_t *f2) {
    return gin_bs_comp(f1->bits, f2->bits);
}

void gin_fmi_serialize_from_buffer(unsigned char *buf, uint64_t buf_size, gin_fmi_t **fmi_ret) {
    gin_bs_t *bs;
    gin_bs_init_from_buffer_copy(buf, buf_size, &bs);
    word_t ridx = 0;
    uint_t fmi_ridx_start = 0;
    word_t fmi_no_bits;
    gin_bs_read_word(bs, ridx, GIN_FMI_NO_BITS_BIT_LENGTH, &fmi_no_bits);
    ridx += GIN_FMI_NO_BITS_BIT_LENGTH;
    gin_fmi_t *fmi = calloc(1, sizeof(gin_fmi_t));
    fmi->no_bits = (int_t)fmi_no_bits;
    /******************************************************
    * Step 1 - Read the header
    ******************************************************/
    word_t no_chars, no_chars_per_block, isa_sample_rate, alphabet_size;
    gin_bs_read_word(bs, ridx, GIN_FMI_CHAR_COUNT_BIT_LENGTH, &no_chars);
    ridx += GIN_FMI_CHAR_COUNT_BIT_LENGTH;

    gin_bs_read_word(bs, ridx, GIN_FMI_RNK_SAMPLE_RATE_BIT_LENGTH, &no_chars_per_block);
    ridx += GIN_FMI_RNK_SAMPLE_RATE_BIT_LENGTH;

    gin_bs_read_word(bs, ridx, GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH, &isa_sample_rate);
    ridx += GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;

    gin_bs_read_word(bs, ridx, GIN_FMI_ALPHABET_SIZE_BIT_LENGTH, &alphabet_size);
    ridx += GIN_FMI_ALPHABET_SIZE_BIT_LENGTH;

    fmi->no_chars = (int_t)no_chars;
    fmi->no_chars_per_block = (int_t)no_chars_per_block;
    fmi->isa_sample_rate = (int_t)isa_sample_rate;
    fmi->alphabet_size = (int_t)alphabet_size;
    fmi->no_bits_per_char = gin_ceil_log2((int_t)alphabet_size);
    /******************************************************
    * Step 2 - Read the alphabet
    ******************************************************/
    fmi->alphabet = calloc(fmi->alphabet_size, sizeof(int_t));
    fmi->e2c = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    fmi->c2e = calloc(GIN_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    for(int_t i = 0; i < fmi->alphabet_size; i++) {
        word_t alphabet_char, encoding;
        gin_bs_read_word(bs, ridx, GIN_FMI_ALPHABET_ENTRY_BIT_LENGTH, &alphabet_char);
        ridx+=GIN_FMI_ALPHABET_ENTRY_BIT_LENGTH;
        fmi->alphabet[i] = (int_t)alphabet_char;

        gin_bs_read_word(bs, ridx, GIN_FMI_ALPHABET_ENCODING_BIT_LENGTH, &encoding);
        ridx+=GIN_FMI_ALPHABET_ENCODING_BIT_LENGTH;

        fmi->c2e[alphabet_char] = (int_t)encoding;
        fmi->e2c[encoding] = (int_t)alphabet_char;
    }
    /****************************
    * Align to word boundary
    ****************************/
    ridx=fmi_ridx_start+((1+(((ridx-fmi_ridx_start)-1)>>WORD_LOG_BITS))<<WORD_LOG_BITS);
    /******************************************************
    * Step 3 - Copy the bits field, then set pointers to
    * suffix array entries and occupancy bitvector
    ******************************************************/
    gin_bs_t *fmi_bits;
    gin_bs_init(&fmi_bits);
    gin_bs_fit(fmi_bits, fmi_no_bits);
    fmi->sa_start_offset = (int_t)(ridx) - (int_t)fmi_ridx_start;
    ridx += (1+fmi->no_chars/fmi->isa_sample_rate)*GIN_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;
    fmi->sa_bv_start_offset = (int_t)ridx - (int_t)fmi_ridx_start;
    int_t sa_bv_occ_size = (1+(fmi->no_chars-1)/ GIN_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH) << GIN_FMI_SA_OCC_BV_LOG_BLOCK_SIZE;
    ridx += sa_bv_occ_size;
    /****************************
    * Align to word boundary
    ****************************/
    ridx=fmi_ridx_start+((1+(((ridx-fmi_ridx_start)-1)>>WORD_LOG_BITS))<<WORD_LOG_BITS);
    /******************************************************
    * Step 4 - Copy the actual bitvector and the caches
    ******************************************************/
    fmi->bv_start_offset = (int_t)ridx - (int_t)fmi_ridx_start;
    for(int_t i = 0; i < fmi_bits->cap_in_words; i++) {
        word_t read_bits;
        uint_t offset = i * WORD_NUM_BITS;
        gin_bs_read_word(bs, fmi_ridx_start + offset, WORD_NUM_BITS, &read_bits);
        gin_bs_write_word(fmi_bits, offset, read_bits, WORD_NUM_BITS);
    }
    gin_bs_fit(fmi_bits, fmi_no_bits);
    fmi->bits = fmi_bits;
    fmi->no_bits = (int_t)fmi_no_bits;
    /******************************************************
    * Step 5 - Compute the cumulative sum from last block
    ******************************************************/
    fmi->char_counts = calloc(fmi->alphabet_size, sizeof(count_t));
    count_t cum = 0;
    uint_t rank_cache_idx = fmi->no_bits - fmi->alphabet_size * GIN_FMI_CHAR_COUNT_BIT_LENGTH;
    for(int_t i = 0; i < fmi->alphabet_size; i++) {
        word_t count;
        fmi->char_counts[i] = cum;
        gin_bs_read_word(fmi->bits, rank_cache_idx + i * GIN_FMI_CHAR_COUNT_BIT_LENGTH, GIN_FMI_CHAR_COUNT_BIT_LENGTH, &count);
        cum += (count_t)count;
    }
    /******************************************************
    * Step 6 - Finalize reading the FMI
    ******************************************************/
    gin_bs_free(bs);
    *fmi_ret = fmi;
}

void gin_fmi_serialize_to_buffer(gin_fmi_t *fmi, unsigned char **buf_ret, uint64_t *buf_size_re) {
    if(!fmi->no_bits) {
        *buf_ret = NULL;
        *buf_size_re = 0;
        return;
    }
    gin_bs_fit(fmi->bits, fmi->no_bits);
    gin_bs_t *copy = gin_bs_copy(fmi->bits);
    *buf_ret = (unsigned char*)copy->words;
    *buf_size_re = 1 + ((fmi->no_bits - 1) >> 3);
    gin_bs_free_disown(copy);
}

void gin_fmi_bwt(gin_fmi_t *fmi, uint64_t *buf, uint64_t start, uint64_t end) {
    for(uint64_t i = start; i <= end; i++) {
        buf[i] = gin_fmi_get(fmi, (int64_t)i);
    }
}