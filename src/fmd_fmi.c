/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_fmi.c is part of fmd
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
#include "fmd_fmi.h"
#include "string.h"

fmd_fmi_qr_t *fmd_fmi_qr_init(int_t lo, int_t hi, int_t pos, fmd_string_t *pattern) {
    fmd_fmi_qr_t *qr = calloc(1, sizeof(fmd_fmi_qr_t));
    qr->lo = lo;
    qr->hi = hi;
    qr->pos = pos;
    qr->pattern = pattern;
    return qr;
}
void fmd_fmi_qr_free(fmd_fmi_qr_t *qr) {
    if(qr) free(qr);
}

void fmd_fmi_init(fmd_fmi_t **fmi,
                  fmd_string_t *string,
                  int_t rank_sample_rate,
                  int_t isa_sample_rate) {
    int64_t *sa = calloc(string->size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)string->seq, (saidx64_t*)sa, string->size+1);
    fmd_fmi_init_with_sa(fmi,string,sa,rank_sample_rate,isa_sample_rate);
    free(sa);
}

void fmd_fmi_init_with_sa(fmd_fmi_t **fmi,
                          fmd_string_t *string,
                          int_t *sa,
                          int_t rank_sample_rate,
                          int_t isa_sample_rate) {
    fmd_fmi_t *f = calloc(1, sizeof(fmd_fmi_t));
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
    char occ[FMD_FMI_MAX_ALPHABET_SIZE];
    memset(occ, 0, FMD_FMI_MAX_ALPHABET_SIZE);
    for(int_t i = 0; i < string->size; i++) {
        char_t c = string->seq[i];
        occ[c] = 1;
    }
    /******************************************************
    * Step 1b - Fill encoding and decoding tables
    ******************************************************/
    int_t alphabet_size = 0;
    for(int_t i = 0; i < FMD_FMI_MAX_ALPHABET_SIZE; i++) {
        alphabet_size += occ[i];
    }
    f->alphabet_size = alphabet_size;
    f->no_bits_per_char = fmd_ceil_log2(f->alphabet_size);
    f->alphabet = calloc(alphabet_size, sizeof(int_t));
    f->e2c = calloc(FMD_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    f->c2e = calloc(FMD_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    int_t t = 0;
    for(int_t i = 0; i < FMD_FMI_MAX_ALPHABET_SIZE; i++) {
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
    fmd_tree_t *sa_tree;
    fmd_tree_init(&sa_tree, &prm_fstruct, &prm_fstruct); // tree sort
    if(!sa_tree) {
        free(f->alphabet);
        free(f->c2e);
        free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    fmd_string_t *bwt;
    fmd_string_init(&bwt, (int_t)f->no_chars);
    for(int_t i = 0; i < (int_t)f->no_chars; i++) {
        int64_t sa_val = sa[i];
        bwt->seq[i] = sa_val ? string->seq[sa_val-1] : FMD_STRING_TERMINATOR;
        if((sa_val % f->isa_sample_rate) == 0) {
            fmd_tree_insert(sa_tree, (void*)i, (void*)sa_val);
        }
    }
    bwt->size = string->size+1;
    /**************************************************************************
    * Step 3 - Encode everything into the bitvector (blocks, ranks, sa, abc)
    *************************************************************************/
    fmd_bs_t *bits;
    fmd_bs_init(&bits);
    if(!bits) {
        fmd_string_free(bwt);
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
    int_t widx = 0;
    fmd_bs_write_word(bits, widx, f->no_chars, FMD_FMI_CHAR_COUNT_BIT_LENGTH);
    widx+=FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    fmd_bs_write_word(bits, widx, f->no_chars_per_block, FMD_FMI_RNK_SAMPLE_RATE_BIT_LENGTH);
    widx+=FMD_FMI_RNK_SAMPLE_RATE_BIT_LENGTH;
    fmd_bs_write_word(bits, widx, f->isa_sample_rate, FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH);
    widx+=FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH;
    fmd_bs_write_word(bits, widx, f->alphabet_size, FMD_FMI_ALPHABET_SIZE_BIT_LENGTH);
    widx+=FMD_FMI_ALPHABET_SIZE_BIT_LENGTH;
    /******************************************************
    * Step 3b - Write the alphabet
    ******************************************************/
    for(int_t i = 0; i < f->alphabet_size; i++) {
        fmd_bs_write_word(bits, widx, (word_t)f->alphabet[i], FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH);
        widx+=FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH;
        word_t encoding = f->c2e[f->alphabet[i]];
        fmd_bs_write_word(bits, widx, encoding, FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH);
        widx+=FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH;
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
    widx += FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH*(1+f->no_chars/f->isa_sample_rate);
    f->sa_bv_start_offset = widx;
    int_t sa_bv_no_blocks = 1 + (f->no_chars - 1) / FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH;
    // reserve space for the bitvector with a dud write
    fmd_bs_write_word(bits, f->sa_bv_start_offset + (sa_bv_no_blocks << FMD_FMI_SA_OCC_BV_LOG_BLOCK_SIZE), 0, 1);
    fmd_fmi_init_isa_write_p_t write_p;
    write_p.bits = bits;
    write_p.sa_start_offset = f->sa_start_offset;
    write_p.sa_bv_start_offset = f->sa_bv_start_offset;
    write_p.isa_sampling_rate = f->isa_sample_rate;
    write_p.no_traversed = 0;
    fmd_tree_inorder(sa_tree, &write_p, fmd_fmi_init_isa_write);
    fmd_tree_free(sa_tree);

    // write the popcounts
    int_t sa_ridx = f->sa_bv_start_offset;
    uint_t sa_bv_cur_popcnt = 0;
    int_t no_words_per_block = FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH >> WORD_LOG_BITS;
    for(int_t i = 0; i < sa_bv_no_blocks; i++) {
        fmd_bs_write_word(bits, sa_ridx, sa_bv_cur_popcnt, FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH);
        sa_ridx += FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH;
        for(int_t j = 0; j < no_words_per_block; j++) {
            word_t payload = 0;
            fmd_bs_read_word(bits, sa_ridx, WORD_NUM_BITS, &payload);
            sa_ridx += WORD_NUM_BITS;
            sa_bv_cur_popcnt += fmd_popcount64(payload);
        }
    }
    int_t sa_bv_occ_size = sa_bv_no_blocks << FMD_FMI_SA_OCC_BV_LOG_BLOCK_SIZE;
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
        fmd_bs_write_word(bits, widx, block_count[c], FMD_FMI_CHAR_COUNT_BIT_LENGTH);
        widx += FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    }
    for(int_t i = 0; i < no_blocks; i++) {
        int_t s = i * (int_t)f->no_chars_per_block;
        int_t e = MIN2((i + 1) * (int_t)f->no_chars_per_block,bwt->size);
        for(int_t j = s; j < e; j++) {
            int_t enc = f->c2e[bwt->seq[j]];
            ++block_count[enc];
            fmd_bs_write_word(bits, widx, enc, f->no_bits_per_char);
            widx += f->no_bits_per_char;
        }
        for(int_t c = 0; c < f->alphabet_size; c++) {
            fmd_bs_write_word(bits, widx, block_count[c], FMD_FMI_CHAR_COUNT_BIT_LENGTH);
            widx += FMD_FMI_CHAR_COUNT_BIT_LENGTH;
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
    fmd_bs_fit(bits, widx);
    f->bits = bits;
    f->no_bits = widx;
    fmd_string_free(bwt);
    *fmi = f;
}

count_t fmd_fmi_query_count(fmd_fmi_t *fmi, fmd_string_t *string) {
    fmd_fmi_qr_t query;
    query.lo = 0;
    query.hi = fmi->no_chars+1;
    query.pos = string->size-1;
    query.pattern = string;
    while(query.pos > -1) {
        bool non_empty = fmd_fmi_advance_query(fmi, &query);
        if(!non_empty) return 0;
    }
    return query.hi - query.lo;
}

fmd_vector_t *fmd_fmi_query_locate(fmd_fmi_t *fmi, fmd_string_t *string) {
    fmd_fmi_qr_t query;
    query.lo = 0;
    query.hi = fmi->no_chars+1;
    query.pos = string->size-1;
    query.pattern = string;
    while(query.pos > -1) {
        bool non_empty = fmd_fmi_advance_query(fmi, &query);
        if(!non_empty) break;
    }
    return fmd_fmi_sa(fmi, &query);
}

void fmd_fmi_init_isa_write(void *key, void *value, void *params) {
    // key: rank, value: offset
    int_t rank = (int_t)key;
    int_t offset = (int_t)value;
    fmd_fmi_init_isa_write_p_t *p = (fmd_fmi_init_isa_write_p_t*)params;
    uint_t offset_widx = p->no_traversed * FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH + p->sa_start_offset;
    fmd_bs_write_word(p->bits, offset_widx, offset, FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH);

    int_t div = rank / FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH; // div is the index of the block
    int_t rem = rank % FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH;
    uint_t sa_occ_widx = p->sa_bv_start_offset + FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH + (div << FMD_FMI_SA_OCC_BV_LOG_BLOCK_SIZE) + rem;
    fmd_bs_write_word(p->bits, sa_occ_widx, (word_t)1, 1);
    p->no_traversed++;
}

bool fmd_fmi_advance_query(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    encoding = fmi->c2e[qr->pattern->seq[qr->pos]];
    count_t rank_lo_m_1 = qr->lo ? fmd_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? fmd_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    qr->lo = (int_t)(base + rank_lo_m_1);
    qr->hi = (int_t)(base + rank_hi_m_1);
    --qr->pos;
    return true;
}

bool fmd_fmi_query_precedence_range(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr, char_t c, int_t *lo, int_t *hi) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    encoding = fmi->c2e[qr->pattern->seq[qr->pos]];
    count_t rank_lo_m_1 = qr->lo ? fmd_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? fmd_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    *lo = (int_t)(base + rank_lo_m_1);
    *hi = (int_t)(base + rank_hi_m_1);
    return true;
}

count_t fmd_fmi_rank(fmd_fmi_t *fmi, word_t enc, int_t pos) {
    // todo: read from short side
    uint_t block_idx = pos / fmi->no_chars_per_block;
    uint_t block_size = (fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    uint_t rank_cache_idx = fmi->bv_start_offset + block_idx * block_size;
    uint_t chars_idx = rank_cache_idx + fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    uint_t char_offset_into_block = pos % fmi->no_chars_per_block;
    word_t rank = 0;
    fmd_bs_read_word(fmi->bits, rank_cache_idx + enc * FMD_FMI_CHAR_COUNT_BIT_LENGTH, FMD_FMI_CHAR_COUNT_BIT_LENGTH, &rank);
    for(int_t i = 0; i <= char_offset_into_block; i++) {
        word_t read_enc = 0;
        fmd_bs_read_word(fmi->bits, chars_idx, fmi->no_bits_per_char, &read_enc);
        rank += enc == read_enc;
        chars_idx += fmi->no_bits_per_char;
    }
    return rank;
}

fmd_vector_t *fmd_fmi_sa(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr) {
    fmd_vector_t *rval;
    fmd_vector_init(&rval, (int_t)qr->hi - (int_t)qr->lo, &prm_fstruct);
    if(!rval) return NULL;
    #pragma omp parallel for default(none) shared(qr, fmi, rval)
    for(int_t j = (int_t)qr->lo; j < (int_t)qr->hi; j++) {
        count_t count = 0;
        int_t i = j;
        while(1) {
            // compute physical index of the SA entry in the bv
            int_t div = i / (FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH);
            int_t rem = i % (FMD_FMI_SA_OCC_BV_PAYLOAD_BIT_LENGTH);
            int_t sa_idx = fmi->sa_bv_start_offset + FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH + (div << FMD_FMI_SA_OCC_BV_LOG_BLOCK_SIZE) + rem;
            word_t sampled = 0;
            fmd_bs_read_word(fmi->bits, sa_idx, 1, &sampled);
            if(!sampled) {
                // sa sample not found, traverse bwt
                count++;
                // get rank of char at index i of the bwt
                word_t encoding = fmd_fmi_get(fmi, (int_t)i);
                // get rank of the char at the same index
                uint64_t rank = fmd_fmi_rank(fmi, encoding, (int_t)i);
                // compute the next index
                i = (int_t)(fmi->char_counts[encoding] + rank - 1);
            } else {
                // sa entry is present, interpolate and compute entry in next range
                word_t popcnt = 0;
                word_t _ = 0;
                word_t sa_entry = 0;
                // read the popcount of the sa occupancy bitvector at the position
                int_t block_base = fmi->sa_bv_start_offset + (div << FMD_FMI_SA_OCC_BV_LOG_BLOCK_SIZE);
                fmd_bs_read_word(fmi->bits,
                                 block_base,
                                 FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH,
                                 &popcnt);
                int_t no_words_before_pos = (rem+1) >> WORD_LOG_BITS;
                int_t no_slack_bits = (rem+1)&(int_t)WORD_LOG_MASK;// ? rem & (int_t)WORD_LOG_MASK + 1 : 0;
                int_t in_block_popcnt_ridx = block_base + FMD_FMI_SA_OCC_BV_POPCOUNT_BIT_LENGTH;
                for(int_t k = 0; k < no_words_before_pos; k++) {
                    _ = 0;
                    fmd_bs_read_word(fmi->bits,
                                     in_block_popcnt_ridx,
                                     WORD_NUM_BITS,
                                     &_);
                    popcnt += fmd_popcount64(_);
                    in_block_popcnt_ridx += WORD_NUM_BITS;
                }
                _ = 0;
                fmd_bs_read_word(fmi->bits,
                                 in_block_popcnt_ridx,
                                 no_slack_bits,
                                 &_);
                popcnt += fmd_popcount64(_);
                // read the SA entry at index popcnt
                fmd_bs_read_word(fmi->bits,
                                 fmi->sa_start_offset + (popcnt-1)*FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH,
                                 FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH,
                                 &sa_entry);

                rval->data[j - qr->lo] = (void*)(sa_entry + count);
                break;
            }
        }
    }
    rval->size = (int_t)qr->hi - (int_t)qr->lo;
    return rval;
}

word_t fmd_fmi_get(fmd_fmi_t *fmi, int_t pos) {
    // todo: read from short side
    uint_t block_idx = pos / fmi->no_chars_per_block;
    uint_t block_size = (fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    uint_t char_idx = fmi->bv_start_offset + block_idx * block_size + fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + (pos % fmi->no_chars_per_block) * fmi->no_bits_per_char;
    word_t read_enc = 0;
    fmd_bs_read_word(fmi->bits, char_idx, fmi->no_bits_per_char, &read_enc);
    return read_enc;
}

fmd_fmi_t* fmd_fmi_copy(fmd_fmi_t* fmi) {
    fmd_fmi_t *f = calloc(1, sizeof(fmd_fmi_t));
    f->no_chars = fmi->no_chars;
    f->no_chars_per_block = fmi->no_chars_per_block;
    f->isa_sample_rate = fmi->isa_sample_rate;
    f->no_bits_per_char = fmi->no_bits_per_char;
    f->alphabet_size = fmi->alphabet_size;
    f->alphabet = calloc(fmi->alphabet_size, sizeof(int_t));
    memcpy(f->alphabet, fmi->alphabet, sizeof(int_t) * fmi->alphabet_size);
    f->c2e = calloc(FMD_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    memcpy(f->c2e, fmi->c2e, sizeof(int_t) * FMD_FMI_MAX_ALPHABET_SIZE);
    f->e2c = calloc(FMD_FMI_MAX_ALPHABET_SIZE, sizeof(int_t));
    memcpy(f->e2c, fmi->e2c, sizeof(int_t) * FMD_FMI_MAX_ALPHABET_SIZE);
    f->char_counts = calloc(f->alphabet_size, sizeof(count_t));
    memcpy(f->char_counts, fmi->char_counts, sizeof(count_t)*f->alphabet_size);
    f->bv_start_offset = fmi->bv_start_offset;
    f->bits = fmd_bs_copy(fmi->bits);
    f->no_bits = fmi->no_bits;
    return f;
}

void fmd_fmi_free(fmd_fmi_t *fmi) {
    if(fmi) {
        free(fmi->alphabet);
        free(fmi->c2e);
        free(fmi->e2c);
        free(fmi->char_counts);
        fmd_bs_free(fmi->bits);
        free(fmi);
    }
}

void fmd_fmi_free_disown(fmd_fmi_t *fmi) {
    if(fmi) {
        free(fmi->alphabet);
        free(fmi->c2e);
        free(fmi->e2c);
        free(fmi->char_counts);
        fmd_bs_free_disown(fmi->bits);
        free(fmi);
    }
}

uint_t fmd_fmi_hash(fmd_fmi_t *fmi) {
    return fmd_bs_hash(fmi->bits);
}

int fmd_fmi_comp(fmd_fmi_t *f1, fmd_fmi_t *f2) {
    return fmd_bs_comp(f1->bits, f2->bits);
}