/*
 * fmd: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023, Unsal Ozturk
 *
 * fmd_bitstream.c is part of fmd
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

void fmd_bs_init_from_buffer(unsigned char *buf, size_t buf_size, fmd_bs_t **bs) {
    fmd_bs_t *b = calloc(1, sizeof(fmd_bs_t));
    if(!b || !buf_size) {
        *bs = NULL;
        return;
    }
    int_t no_words = (1 + ((int_t)buf_size - 1) / (int_t)sizeof(word_t));
    b->words = realloc(buf, sizeof(word_t) * no_words);
    memset((unsigned char*)b->words + buf_size, 0, no_words * sizeof(word_t) - buf_size);
    if(!b->words) {
        free(b);
        *bs = NULL;
        return;
    }
    b->cap_in_words = no_words;
    *bs = b;
}

void fmd_bs_read_word(fmd_bs_t *bs, uint_t bit_idx, uint_t read_size_in_bits, word_t *read_val) {
    if(!read_size_in_bits) return;
    uint_t word_idx = bit_idx >> WORD_LOG_BITS;
    uint_t bit_s = (bit_idx & WORD_LOG_MASK);
    uint_t bit_e = (read_size_in_bits + bit_idx - 1) & WORD_LOG_MASK;
    if (bit_s == 0) {
        *read_val = bs->words[word_idx] & mask_shift_right_64[WORD_NUM_BITS - read_size_in_bits];
    } else if (bit_e >= bit_s) {
        *read_val = (bs->words[word_idx] & (mask_shift_left_64[bit_s] & mask_shift_right_64[WORD_NUM_BITS - 1 - bit_e])) >> bit_s;
    } else {
        *read_val = ((bs->words[word_idx] & mask_shift_left_64[bit_s]) >> (bit_s)) |
                    ((bs->words[word_idx + 1] & mask_shift_right_64[WORD_NUM_BITS - 1 - bit_e]) << (WORD_NUM_BITS - bit_s));
    }
}

void fmd_bs_write_word(fmd_bs_t *bs, uint_t bit_idx, word_t write_val, uint_t write_size_in_bits) {
    if(write_size_in_bits == 0) return;
    while (bit_idx + write_size_in_bits > (bs->cap_in_words << WORD_LOG_BITS)) {
        bs->words = realloc(bs->words, sizeof(word_t) * bs->cap_in_words * 2);
        memset(bs->words + bs->cap_in_words, 0, bs->cap_in_words * sizeof(word_t));
        bs->cap_in_words *= 2;
    }
    uint_t word_idx = bit_idx >> WORD_LOG_BITS;
    uint_t bit_s = bit_idx & WORD_LOG_MASK;
    uint_t bit_e = (bit_idx + write_size_in_bits - 1) & WORD_LOG_MASK;
    if (bit_s == 0) { // word aligned write, write no_bits_to_write many bits from data
        bs->words[word_idx] &= mask_shift_left_64[write_size_in_bits];
        bs->words[word_idx] |= write_val & mask_shift_right_64[WORD_NUM_BITS - write_size_in_bits];
    } else if (bit_e >= bit_s) { // write into a single word, clear the middle bits
        word_t mask = mask_shift_left_64[bit_s] & mask_shift_right_64[WORD_NUM_BITS - 1 - bit_e];
        bs->words[word_idx] &= ~mask;
        bs->words[word_idx] |= mask & (write_val << (bit_s));
    } else { // write to two adjacent words
        bs->words[word_idx] &= ~mask_shift_left_64[bit_s];
        bs->words[word_idx++] |= (write_val << bit_s);
        bs->words[word_idx] &= mask_shift_left_64[++bit_e];
        bs->words[word_idx] |= (write_val >> (WORD_NUM_BITS - bit_s)) & ~mask_shift_left_64[bit_e];
    }
}

void fmd_bs_fit(fmd_bs_t *bs, uint_t bit_idx) {
    uint_t final_cap = 1 + (uint_t)((int_t)bit_idx - 1) / WORD_NUM_BITS;
    uint_t prev_cap = bs->cap_in_words;
    if(final_cap == prev_cap) return;
    bs->words = realloc(bs->words, final_cap * sizeof(word_t));
    if(final_cap > prev_cap) {
        memset(bs->words + prev_cap, 0, (final_cap - prev_cap) * sizeof(word_t));
    }
    bs->cap_in_words = final_cap;
}

void fmd_bs_detach(fmd_bs_t *bs, word_t **words, uint_t *no_words) {
    *words = bs->words;
    *no_words = bs->cap_in_words;
    bs->words = calloc(BS_INIT_SIZE_WORDS, sizeof(word_t));
    if(!bs->words) {
        free(bs);
        bs = NULL;
        return;
    }
    bs->cap_in_words = BS_INIT_SIZE_WORDS;
}

void fmd_bs_free(fmd_bs_t *bs) {
    if(bs) {
        free(bs->words);
        free(bs);
    }
}

void fmd_bs_free_disown(fmd_bs_t *bs) {
    if(bs) {
        free(bs);
    }
}

fmd_bs_t *fmd_bs_copy(fmd_bs_t *bs) {
    fmd_bs_t *cpy;
    fmd_bs_init_reserve(&cpy, bs->cap_in_words);
    if(!cpy) return NULL;
    memcpy(cpy->words, bs->words, sizeof(word_t)*bs->cap_in_words);
    return cpy;
}

int fmd_bs_comp(fmd_bs_t *bs1, fmd_bs_t *bs2) {
    return bs1->cap_in_words == bs2->cap_in_words ? memcmp(bs1->words, bs2->words, bs1->cap_in_words * sizeof(word_t)) : -1;
}

uint_t fmd_bs_hash(fmd_bs_t *bs) {
    // fnv1a_64
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;

    const uint8_t* data = (const uint8_t*)bs->words;
    for (size_t i = 0; i < bs->cap_in_words; ++i) {
        hash ^= (uint64_t)data[i];
        hash *= prime;
    }
    return hash;
}