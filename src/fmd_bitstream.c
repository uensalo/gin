#include "../include/fmd_bitstream.h"

void fmd_bs_init(fmd_bs_t **bs) {
    fmd_bs_t *b = calloc(1, sizeof(fmd_bs_t));
    if(!b) {
        *bs = NULL;
        return;
    }
    b->words = calloc(BS_INIT_SIZE_WORDS, sizeof(word_t));
    if(!b->words) {
        free(b);
        *bs = NULL;
        return;
    }
    b->cap_in_words = BS_INIT_SIZE_WORDS;
    *bs = b;
}

void fmd_bs_init_reserve(fmd_bs_t **bs, uint64_t no_words) {
    fmd_bs_t *b = calloc(1, sizeof(fmd_bs_t));
    if(!b) {
        *bs = NULL;
        return;
    }
    b->words = calloc(no_words, sizeof(word_t));
    if(!b->words) {
        free(b);
        *bs = NULL;
        return;
    }
    b->cap_in_words = no_words;
    *bs = b;
}

void fmd_bs_init_from_buffer(char *buf, size_t buf_size, fmd_bs_t **bs) {
    fmd_bs_t *b = calloc(1, sizeof(fmd_bs_t));
    if(!b) {
        *bs = NULL;
        return;
    }
    b->words = calloc(1, sizeof(word_t) * ((buf_size >> WORD_LOG_BITS) + 1));
    if(!b->words) {
        free(b);
        *bs = NULL;
        return;
    }
    memcpy((char*)b->words, buf, buf_size);
    *bs = b;
}

void fmd_bs_init_read_word(fmd_bs_t *bs, upos_t bit_idx, upos_t read_size_in_bits, word_t *read_val) {
    if(!read_size_in_bits) return;
    upos_t word_idx = bit_idx >> WORD_LOG_BITS;
    upos_t bit_s = (bit_idx & WORD_LOG_MASK);
    upos_t bit_e = (read_size_in_bits + bit_idx - 1) & WORD_LOG_MASK;
    if (bit_s == 0) {
        *read_val = bs->words[word_idx] & mask_shift_right_64[64 - read_size_in_bits];
    } else if (bit_e >= bit_s) {
        *read_val = (bs->words[word_idx] & (mask_shift_left_64[bit_s] & mask_shift_right_64[63 - bit_e])) >> bit_s;
    } else {
        *read_val = ((bs->words[word_idx] & mask_shift_left_64[bit_s]) >> (bit_s)) |
                    ((bs->words[word_idx + 1] & mask_shift_right_64[63 - bit_e]) << (64 - bit_s));
    }
}

void fmd_bs_init_write_word(fmd_bs_t *bs, upos_t bit_idx, word_t write_val, upos_t write_size_in_bits) {
    if(write_size_in_bits == 0) return;
    while (bit_idx + write_size_in_bits > (bs->cap_in_words << WORD_LOG_BITS)) {
        bs->words = realloc(bs->words, sizeof(word_t) * bs->cap_in_words * 2);
        memset(bs->words + bs->cap_in_words, 0, bs->cap_in_words * sizeof(word_t));
        bs->cap_in_words *= 2;
    }
    upos_t word_idx = bit_idx >> WORD_LOG_BITS;
    upos_t bit_s = bit_idx & WORD_LOG_MASK;
    upos_t bit_e = (bit_idx + write_size_in_bits - 1) & WORD_LOG_MASK;
    if (bit_s == 0) { // word aligned write, write no_bits_to_write many bits from data
        bs->words[word_idx] &= mask_shift_left_64[write_size_in_bits];
        bs->words[word_idx] |= write_val & mask_shift_right_64[64 - write_size_in_bits];
    } else if (bit_e >= bit_s) { // write into a single word, clear the middle bits
        word_t mask = mask_shift_left_64[bit_s] & mask_shift_right_64[63 - bit_e];
        bs->words[word_idx] &= ~mask;
        bs->words[word_idx] |= mask & (write_val << (bit_s));
    } else { // write to two adjacent words
        bs->words[word_idx] &= ~mask_shift_left_64[bit_s];
        bs->words[word_idx++] |= (write_val << bit_s);
        bs->words[word_idx] &= mask_shift_left_64[++bit_e];
        bs->words[word_idx] |= (write_val >> (64 - bit_s)) & ~mask_shift_left_64[bit_e];
    }
}

void fmd_bs_free(fmd_bs_t *bs) {
    if(bs) {
        free(bs->words);
        free(bs);
    }
}