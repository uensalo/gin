#include "fmd_fmi.h"
#include "string.h"

void fmd_fmi_init(fmd_fmi_t **fmi,
                  fmd_string_t *string,
                  pos_t rank_sample_rate,
                  pos_t isa_sample_rate) {
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
    * Step 1a - Initialize tables and alphabet
    ******************************************************/
    fmd_vector_init(&f->alphabet, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    if(!f->alphabet) {
        free(f);
        *fmi = NULL;
        return;
    }
    // init tables
    fmd_table_init(&f->c2e, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    if(!f->c2e) {
        fmd_vector_free(f->alphabet);
        free(f);
        *fmi = NULL;
        return;
    }
    fmd_table_init(&f->e2c, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    if(!f->e2c) {
        fmd_vector_free(f->alphabet);
        fmd_table_free(f->c2e);
        free(f);
        *fmi = NULL;
        return;
    }
    /******************************************************
    * Step 1b - Insert unique characters into a tree
    ******************************************************/
    fmd_tree_t *char_set;
    fmd_tree_init(&char_set, &prm_fstruct, &prm_fstruct);
    fmd_tree_insert(char_set, (void*)FMD_STRING_TERMINATOR, (void*)1);
    if(!char_set) {
        fmd_vector_free(f->alphabet);
        fmd_table_free(f->c2e);
        fmd_table_free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    for(pos_t i = 0; i < string->size; i++) {
        char_t c = string->seq[i];
        fmd_tree_insert(char_set, (void*)c, 0);
    }
    /******************************************************
    * Step 1c - Flatten tree into a sorted list
    ******************************************************/
    fmd_fmi_init_charset_flatten_p_t trav_params;
    trav_params.chars = f->alphabet;
    trav_params.e2c = f->e2c;
    trav_params.c2e = f->c2e;
    fmd_tree_inorder(char_set, &trav_params, fmd_fmi_init_charset_flatten);
    fmd_tree_free(char_set);
    f->alphabet_size = f->alphabet->size;
    f->no_bits_per_char = fmd_ceil_log2(f->alphabet_size);
    /**************************************************************************
    * Step 2 - Compute the suffix array and the BWT of the input and sample
    *************************************************************************/
    fmd_table_init(&f->isa, FMD_HT_INIT_SIZE, &prm_fstruct, &prm_fstruct);
    if(!f->isa) {
        fmd_vector_free(f->alphabet);
        fmd_table_free(f->c2e);
        fmd_table_free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    int64_t *sa = calloc(string->size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)string->seq, (saidx64_t*)sa, string->size+1);
    fmd_string_t *bwt;
    fmd_string_init(&bwt, (pos_t)f->no_chars);
    for(pos_t i = 0; i < (pos_t)f->no_chars; i++) {
        int64_t sa_val = sa[i];
        bwt->seq[i] = sa_val ? string->seq[sa_val-1] : FMD_STRING_TERMINATOR;
        if((sa_val % f->isa_sample_rate) == 0) {
            fmd_table_insert(f->isa, (void*)i, (void*)sa_val);
        }
    }
    bwt->size = string->size+1;
    free(sa);
    /**************************************************************************
    * Step 3 - Encode everything into the bitvector (blocks, ranks, sa, abc)
    *************************************************************************/
    fmd_bs_t *bits;
    fmd_bs_init(&bits);
    if(!bits) {
        fmd_string_free(bwt);
        fmd_vector_free(f->alphabet);
        fmd_table_free(f->c2e);
        fmd_table_free(f->e2c);
        free(f);
        *fmi = NULL;
        return;
    }
    /******************************************************
    * Step 3a - Write the header
    ******************************************************/
    pos_t widx = 0;
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
    for(pos_t i = 0; i < f->alphabet->size; i++) {
        fmd_bs_write_word(bits, widx, (word_t)f->alphabet->data[i], FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH);
        widx+=FMD_FMI_ALPHABET_ENTRY_BIT_LENGTH;
        word_t encoding;
        fmd_table_lookup(f->c2e, f->alphabet->data[i], &encoding);
        fmd_bs_write_word(bits, widx, encoding, FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH);
        widx+=FMD_FMI_ALPHABET_ENCODING_BIT_LENGTH;
    }
    /******************************************************
    * Step 3c - Write suffix array samples
    ******************************************************/
    fmd_fmi_init_isa_write_p_t write_p;
    write_p.bits = bits;
    write_p.base_bit_idx = widx;
    write_p.isa_sampling_rate = f->isa_sample_rate;
    fmd_table_traverse(f->isa, &write_p, fmd_fmi_init_isa_write);
    widx+=FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH*(f->no_chars/f->isa_sample_rate);
    /******************************************************
    * Step 3d - Write rank caches and payloads
    ******************************************************/
    f->bv_start_offset = widx;
    pos_t no_blocks = 1 + (bwt->size - 1) / f->no_chars_per_block;
    count_t *block_count = calloc(f->alphabet_size, sizeof(count_t));
    for(pos_t c = 0; c < f->alphabet_size; c++) {
        fmd_bs_write_word(bits, widx, block_count[c], FMD_FMI_CHAR_COUNT_BIT_LENGTH);
        widx += FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    }
    for(pos_t i = 0; i < no_blocks; i++) {
        pos_t s = i * (pos_t)f->no_chars_per_block;
        pos_t e = MIN2((i + 1) * (pos_t)f->no_chars_per_block,bwt->size);
        for(pos_t j = s; j < e; j++) {
            pos_t enc;
            bool found = fmd_table_lookup(f->c2e, (void*)bwt->seq[j], (void*)&enc);
            ++block_count[enc];
            fmd_bs_write_word(bits, widx, enc, f->no_bits_per_char);
            widx += f->no_bits_per_char;
        }
        for(pos_t c = 0; c < f->alphabet_size; c++) {
            fmd_bs_write_word(bits, widx, block_count[c], FMD_FMI_CHAR_COUNT_BIT_LENGTH);
            widx += FMD_FMI_CHAR_COUNT_BIT_LENGTH;
        }
    }
    /******************************************************
    * Step 3e - Compute a cumulative sum of char counts
    ******************************************************/
    count_t cum = 0;
    f->char_counts = calloc(f->alphabet_size, sizeof(count_t));
    for(pos_t i = 0; i < f->alphabet_size; i++) {
        f->char_counts[i] = cum;
        cum += block_count[i];
    }
    free(block_count);
    f->bits = bits;
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

// helpers
void fmd_fmi_init_charset_flatten(void *key, void *value, void *params) {
    fmd_fmi_init_charset_flatten_p_t *p = (fmd_fmi_init_charset_flatten_p_t*)params;
    word_t encoding = p->chars->size;
    fmd_table_insert(p->c2e, key, (void*)encoding);
    fmd_table_insert(p->e2c, (void*)encoding, key);
    fmd_vector_append(p->chars, key);
}

void fmd_fmi_init_isa_write(void *key, void *value, void *params) {
    fmd_fmi_init_isa_write_p_t *p = (fmd_fmi_init_isa_write_p_t*)params;
    upos_t widx = (upos_t) value / (upos_t) p->isa_sampling_rate;
    upos_t bit_offset = widx * FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH + p->base_bit_idx;
    fmd_bs_write_word(p->bits, bit_offset, (word_t)key, FMD_FMI_ISA_SAMPLE_RATE_BIT_LENGTH);
}

bool fmd_fmi_advance_query(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr) {
    // traverse the LF-mapping
    // compute the rank of the symbol for lo-1 and hi-1
    word_t encoding;
    bool found = fmd_table_lookup(fmi->c2e, (void*)qr->pattern->seq[qr->pos], &encoding);
    if(!found) {
        qr->lo = 0;
        qr->hi = 0;
        //fprintf(stderr,"[fmd_fmi_advance_query]: encoding not found in dictionary, query is NIL\n");
        return false;
    }
    count_t rank_lo_m_1 = qr->lo ? fmd_fmi_rank(fmi,encoding, qr->lo-1) : 0ull;
    count_t rank_hi_m_1 = qr->hi ? fmd_fmi_rank(fmi,encoding, qr->hi-1) : 0ull;
    uint64_t base = fmi->char_counts[encoding];
    qr->lo = (pos_t)(base + rank_lo_m_1);
    qr->hi = (pos_t)(base + rank_hi_m_1);
    --qr->pos;
    return true;
}

count_t fmd_fmi_rank(fmd_fmi_t *fmi, word_t enc, pos_t pos) {
    // todo: read from short side
    upos_t block_idx = pos / fmi->no_chars_per_block;
    upos_t block_size = (fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    upos_t rank_cache_idx = fmi->bv_start_offset + block_idx * block_size;
    upos_t chars_idx = rank_cache_idx + fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH;
    upos_t char_offset_into_block = pos % fmi->no_chars_per_block;
    word_t rank = 0;
    fmd_bs_read_word(fmi->bits, rank_cache_idx + enc * FMD_FMI_CHAR_COUNT_BIT_LENGTH, FMD_FMI_CHAR_COUNT_BIT_LENGTH, &rank);
    for(pos_t i = 0; i <= char_offset_into_block; i++) {
        word_t read_enc = 0;
        fmd_bs_read_word(fmi->bits, chars_idx, fmi->no_bits_per_char, &read_enc);
        rank += enc == read_enc;
        chars_idx += fmi->no_bits_per_char;
    }
    return rank;
}

fmd_vector_t *fmd_fmi_sa(fmd_fmi_t *fmi, fmd_fmi_qr_t *qr) {
    fmd_vector_t *rval;
    fmd_vector_init(&rval, FMD_VECTOR_INIT_SIZE, &prm_fstruct);
    if(!rval) return NULL;
    for(pos_t j = (pos_t)qr->lo; j < (pos_t)qr->hi; j++) {
        count_t count = 0;
        pos_t i = j;
        while(1) {
            count_t sa_entry = -1;
            bool found = fmd_table_lookup(fmi->isa, (void*)i, &sa_entry);
            if(!found) {
                // sa sample not found, traverse bwt
                count++;
                // get rank of char at index i of the bwt
                word_t encoding = fmd_fmi_get(fmi, (pos_t)i);
                // get rank of the char at the same index
                uint64_t rank = fmd_fmi_rank(fmi, encoding, (pos_t)i);
                // compute the next index
                i = (pos_t)(fmi->char_counts[encoding] + rank - 1);
            } else {
                // sa entry is present, interpolate and compute entry in next range
                fmd_vector_append(rval,(void*)(sa_entry + count));
                break;
            }
        }
    }
    return rval;
}

word_t fmd_fmi_get(fmd_fmi_t *fmi, pos_t pos) {
    // todo: read from short side
    upos_t block_idx = pos / fmi->no_chars_per_block;
    upos_t block_size = (fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + fmi->no_bits_per_char * fmi->no_chars_per_block);
    upos_t char_idx = fmi->bv_start_offset + block_idx * block_size + fmi->alphabet_size * FMD_FMI_CHAR_COUNT_BIT_LENGTH + (pos % fmi->no_chars_per_block) * fmi->no_bits_per_char;
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
    f->alphabet = fmd_vector_copy(fmi->alphabet);
    f->c2e = fmd_table_copy(fmi->c2e);
    f->e2c = fmd_table_copy(fmi->e2c);
    f->char_counts = calloc(f->alphabet_size, sizeof(count_t));
    memcpy(f->char_counts, fmi->char_counts, sizeof(count_t)*f->alphabet_size);
    f->isa = fmd_table_copy(fmi->isa);
    f->bv_start_offset = fmi->bv_start_offset;
    f->bits = fmd_bs_copy(fmi->bits);
    return f;
}

void fmd_fmi_free(fmd_fmi_t *fmi) {
    if(fmi) {
        fmd_vector_free(fmi->alphabet);
        fmd_table_free(fmi->c2e);
        fmd_table_free(fmi->e2c);
        free(fmi->char_counts);
        fmd_table_free(fmi->isa);
        fmd_bs_free(fmi->bits);
        free(fmi);
    }
}

upos_t fmd_fmi_hash(fmd_fmi_t *fmi) {
    return fmd_bs_hash(fmi->bits);
}

int fmd_fmi_comp(fmd_fmi_t *f1, fmd_fmi_t *f2) {
    return fmd_bs_comp(f1->bits, f2->bits);
}