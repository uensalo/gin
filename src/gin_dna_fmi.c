/*
 * gin: FM-Index-like graph indexing algorithm toolkit.
 * Copyright (C) 2023-2024, Unsal Ozturk
 *                    2024, Paolo Ribeca
 *
 * gin_dna_fmi.c is part of gin
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
#include "gin_dna_fmi.h"
#include "divsufsort64.h"

gin_dfmi_t* gin_dfmi_build(const char* str, uint64_t size, uint64_t isa_rate) {
    // Run important alignment and size checks, notify user of misalignment if any
    if(sizeof(gin_dfmi_header_t) != 64) {
        fprintf(stderr, "[gin_dna_fmi.c] Warning: gin_dfmi_header_t (size %d) is not 64 bytes in width. This will cause cache line alignment problems.\n", sizeof(gin_dfmi_header_t));
    }
    if(sizeof(gin_dfmi_sb_t) != 128) {
        fprintf(stderr, "[gin_dna_fmi.c] Warning: gin_dfmi_sbc_t (size %d) is not 128 bytes in width. This will cause cache line alignment problems.\n", sizeof(gin_dfmi_sb_t));
    }
    if(sizeof(gin_dfmi_sao_t) != 64) {
        fprintf(stderr, "[gin_dna_fmi.c] Warning: gin_dfmi_sao_t (size %d) is not 64 bytes in width. This will cause cache line alignment problems.\n", sizeof(gin_dfmi_sao_t));
    }
    if(sizeof(gin_uint40_t) != 5) {
        fprintf(stderr, "[gin_dna_fmi.c] Warning: gin_uint40_t (size %d) is not 5 bytes in width. This will cause cache line alignment problems.\n", sizeof(gin_uint40_t));
    }

    gin_dfmi_t *dfmi = calloc(1, sizeof(gin_dfmi_t));
    if(!dfmi) return NULL;
    // Numbers
    uint64_t N = size+1;
    uint64_t n_sb        = DFMI_CEIL(N,192);
    uint64_t n_sa_occ_sb = DFMI_CEIL(N,384);
    uint64_t V = 0; // Number of ( characters
    uint64_t P = 0; // Number of , and . characters
    // Compute the SA and bwt now.
    char *bwt = (char*)malloc(N*sizeof(char));
    int64_t *sa = calloc(size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)str, (saidx64_t*)sa, (saidx64_t)size+1);
    // Compute the BWT
    for(uint64_t i = 0; i < N; i++) {
        if(sa[i]) {
            bwt[i] = str[sa[i] - 1];
        } else {
            dfmi->header.values.str_term_pos = i;
            bwt[i] = (char)'\0';
        }
        V += bwt[i] == '(';
        P += (bwt[i] == ',' || bwt[i] == '.');
    }

    // Count samples manually because there's no way of
    // Knowing if '()' coincides with the sampling period
    // without computing the BWT. The trick is: always store
    // SA entries for '(' and ')'. For permutation chars, don't
    // store anything as they will never be decoded. For the rest,
    // sample at periods.
    uint64_t n_sa_val = 2*V+1;
    for(uint64_t i = 2*V+1+P; i < N; i++) {
        if(!(sa[i] % isa_rate)) {
            n_sa_val++;
        }
    }

    // Number of words needed to store things
    uint64_t nw_header = sizeof(gin_dfmi_header_t) / sizeof(uint64_t);
    uint64_t nw_f = 8;
    uint64_t nw_l = n_sb * sizeof(gin_dfmi_sb_t) / sizeof(uint64_t);
    uint64_t nw_sa_occ = n_sa_occ_sb * sizeof(gin_dfmi_sao_t) / sizeof(uint64_t);
    uint64_t nw_sa = DFMI_CEIL(n_sa_val * sizeof(gin_uint40_t), sizeof(uint64_t));

    // Allocate space for the large buffer containing everything
    uint64_t no_fmi_words =  nw_header + nw_f +  nw_l + nw_sa_occ + nw_sa;
    uint64_t *buf = calloc(no_fmi_words, sizeof(uint64_t));
    uint64_t wptr = 0;

    // populate starts and write simple stuff to the bs
    dfmi->buf = buf;
    dfmi->header.values.fmi_size_in_bytes = no_fmi_words * sizeof(uint64_t);
    dfmi->header.values.isa_rate = isa_rate;
    dfmi->header.values.no_chars = N; // includes terminator
    dfmi->header.values.str_term_pos = -1; // fill in later
    dfmi->header.values.no_sa_values = n_sa_val;
    dfmi->header.values.no_dfmi_x = V;
    wptr += nw_header;


    // Set proper pointers in memory
    dfmi->f      = (dfmi->buf + wptr); wptr += nw_f;
    dfmi->l      = (gin_dfmi_sb_t *)(dfmi->buf + wptr); wptr += nw_l;
    dfmi->sa_occ = (gin_dfmi_sao_t*)(dfmi->buf + wptr); wptr += nw_sa_occ;
    dfmi->sa     = (gin_uint40_t  *)(dfmi->buf + wptr); wptr += nw_sa;
    dfmi->unr    = (dfmi->buf + wptr);

    // Sample the inverse suffix array
    uint64_t k = 0;
    // Always sample the terminator, '(' and ')'.
    for(uint64_t i = 0; i < 2*V+1; i++) {
        DFMI_SAO_SET_BIT(dfmi->sa_occ, i);
        DFMI_UINT40_SET(&dfmi->sa[k], sa[i]);
        k++;
    }
    // Skip the permutation characters and sample the text
    for(uint64_t i =  2*V+1+P; i < N; i++) {
        if(!(sa[i] % isa_rate)) {
            DFMI_SAO_SET_BIT(dfmi->sa_occ, i);
            DFMI_UINT40_SET(&dfmi->sa[k], sa[i]);
            k++;
        }
    }
    // SA is no longer needed, free here
    free(sa);

    // Populate the L bitvector
    for(uint64_t i = 0; i < N; i++) {
        uint64_t enc = gin_dfmi_enc[bwt[i]];
        DFMI_L_SET_BIT(dfmi->l, enc, i); // heavy lifting
    }

    // Tally character counts for L
    uint64_t ct[DFMI_NO_ALL_CHARS];     // Keep global counts
    uint64_t bct[DFMI_NO_ALL_CHARS]; // In superblock-counts
    memset(ct, 0, sizeof(ct));
    memset(bct, 0, sizeof(bct));
    for(uint64_t i = 0; i < N; i++) {
        uint64_t s = i / 192;        // Superblock index
        uint64_t b = (i % 192) / 64; // Block index within the superblock
        if(!(i % 192)) {
            // Reached a new superblock, cache counts now
            for(uint64_t c = 0; c < DFMI_NO_RANKED_CHARS; c++) {
                DFMI_UINT40_SET(&dfmi->l[s].sb.sbc[c], ct[c]);
            }
            // Reset in-superblock counts
            memset(bct, 0, sizeof(bct));
        }
        if(!(i % 64)) {
            // Reached a new block, populate its in-block counts
            for(uint64_t c = 0; c < DFMI_NO_RANKED_CHARS; c++) {
                dfmi->l[s].sb.blocks[b].bc[c] = (uint8_t)bct[c];
            }
        }
        ct[gin_dfmi_enc[bwt[i]]]++;
        bct[gin_dfmi_enc[bwt[i]]]++;
    }

    // Populate F column manually - lexical order of ASCII matters here :(
    // Original lex: \0 ( ) , . A C G N T
    // This is slightly tricky: We won't be storing counts for permutation characters
    // In order, we will store a subset: (,A,C,G,N,T (ASCII alphabetical order)
    dfmi->f[0] = ct[DFMI_$];                                                     // (
    dfmi->f[1] = dfmi->f[0] + ct[DFMI_X] + ct[DFMI_Y] + ct[DFMI_a] + ct[DFMI_b]; // A
    dfmi->f[2] = dfmi->f[1] + ct[DFMI_A];                                        // C
    dfmi->f[3] = dfmi->f[2] + ct[DFMI_C];                                        // G
    dfmi->f[4] = dfmi->f[3] + ct[DFMI_G];                                        // N
    dfmi->f[5] = dfmi->f[4] + ct[DFMI_N];                                        // T

    // Popcount for suffix array occupancies
    uint64_t sao_popcnt = 0;
    uint64_t sao_bct = 0;
    for(uint64_t i = 0; i < n_sa_occ_sb; i++) {
        DFMI_UINT40_SET(&dfmi->sa_occ[i].sa_occ_sbc.popcnt, sao_popcnt);
        for(uint64_t b = 0; b < 5; b++) { // 6: Number of blocks per SAO superblock
            sao_bct += DFMI_POPCOUNT(dfmi->sa_occ[i].sa_occ_sbc.bv[b]);
            dfmi->sa_occ[i].sa_occ_sbc.bc[b] = sao_bct;
        }
        sao_popcnt += DFMI_POPCOUNT(dfmi->sa_occ[i].sa_occ_sbc.bv[5]) + sao_bct; // next block will be cached later
        sao_bct = 0;
    }

    // Regererate the header before returning
    memcpy(dfmi->buf, dfmi->header.buf, sizeof(gin_dfmi_header_t));

    // Free the bwt here
    free(bwt);
    // Done! But at what cost...
    return dfmi;
}


int64_t gin_dfmi_rank(gin_dfmi_t *dfmi, uint64_t pos, char c) {
    // First, find the encoding of the char
    uint8_t enc = gin_dfmi_enc[c];
    // Possible prefetch here -- fetch the entire superblock
    // Then, go to its superblock and retrieve
    // a) superblock count b) block count c) the wavelet
    uint64_t s =  pos / 192;
    uint64_t b  = (pos % 192) / 64;
    uint64_t m  = pos % 64;
    uint64_t sbc = DFMI_UINT40_GET(dfmi->l[s].sb.sbc[enc]);
    uint64_t bc  = (uint64_t)dfmi->l[s].sb.blocks[b].bc[enc];
    uint64_t wav = gin_dfmi_wavelet_enc(enc, dfmi->l[s].sb.blocks[b].bv);
    // Mask the wavelet up until pos, and return the sum of cached counts
    return (int64_t)(sbc + bc + DFMI_POPCOUNT(wav & DFMI_POPCOUNT_MASK(m)));
}

void gin_dfmi_sa(gin_dfmi_t *dfmi, uint64_t *buf, uint64_t start, uint64_t end) {
    for(uint64_t j = start; j <= end; j++) {
        uint64_t count = 0;
        uint64_t i = j;
        while(1) {
            // Check if the SA sample exists by consulting the occupancy bitvector
            uint64_t sampled = 0;
            sampled = DFMI_SAO_GET_BIT(dfmi->sa_occ, i);
            if(!sampled) {
                // SA not sampled, traverse the BWT
                count++;
                // Get encoding of char at index i of L and its rank
                uint8_t enc = DFMI_L_CHAR(dfmi->l,i);
                uint64_t rank = gin_dfmi_rank(dfmi, i, gin_dfmi_dec[enc]);
                // Navigate in the F column to the correct place
                i = (dfmi->f[enc] + rank - 1);
            } else {
                // SA is sampled! Retrieve from the bitvector
                uint64_t sa_idx = DFMI_SAO_RANK(dfmi->sa_occ, i) - 1;
                // Read the SA entry at the correct index
                uint64_t sa_entry = DFMI_UINT40_GET(dfmi->sa[sa_idx]);
                buf[j-start] = sa_entry + count;
                break;
            }
        }
    }
}

void gin_dfmi_bwt(gin_dfmi_t *dfmi, uint64_t *buf, uint64_t start, uint64_t end) {
    for (uint64_t i = start; i <= end; ++i) {
        buf[i - start] = (uint64_t)gin_dfmi_dec[DFMI_L_CHAR(dfmi->l,i)];
    }
}

void gin_dfmi_to_buffer(gin_dfmi_t *dfmi, uint8_t **data, uint64_t *size) {
    // This is the easiest serialisation ever. Seriously.
    *size = dfmi->header.values.fmi_size_in_bytes;
    *data = calloc(1, sizeof(uint8_t) * (*size));
    if(!(*data)) {
        *size = 0;
        return;
    }
    memcpy(*data, dfmi->buf, *size);
}

void* gin_dfmi_from_buffer(uint8_t *data, uint64_t size) {
    if(!size) return NULL;
    // The buffer will be padded to 64 bit words in a future patch
    //if(size % 64) {
    //    fprintf(stderr, "[gin_dna_fmi] Warning: gin_dna_fmi buffer not aligned to 64 bytes!\n");
    //}
    gin_dfmi_t *dfmi = calloc(1,sizeof(gin_dfmi_t));
    if(!dfmi) {
        return NULL;
    }
    // Copy the data
    uint64_t *copy = calloc(1, sizeof(uint8_t) * size);
    if(!copy) {
        free(dfmi);
        return NULL;
    }
    memcpy(copy, data, size * sizeof(uint8_t));
    // Parse the header
    dfmi->buf = (uint64_t*)copy;
    dfmi->header = *(gin_dfmi_header_t*)dfmi->buf;
    // Calculate sizes
    uint64_t N = dfmi->header.values.no_chars;
    uint64_t n_sb        = DFMI_CEIL(N,192);
    uint64_t n_sa_occ_sb = DFMI_CEIL(N,384);
    uint64_t n_sa_val =    DFMI_CEIL(N, dfmi->header.values.no_chars);
    // Number of words needed to store things
    uint64_t nw_header = sizeof(gin_dfmi_header_t) / sizeof(uint64_t);
    uint64_t nw_f = 8;
    uint64_t nw_l = n_sb * sizeof(gin_dfmi_sb_t) / sizeof(uint64_t);
    uint64_t nw_sa_occ = n_sa_occ_sb * sizeof(gin_dfmi_sao_t) / sizeof(uint64_t);
    uint64_t nw_sa = DFMI_CEIL(n_sa_val * sizeof(gin_uint40_t), sizeof(uint64_t));
    uint64_t wptr = nw_header;
    dfmi->f      = (dfmi->buf + wptr); wptr += nw_f;
    dfmi->l      = (gin_dfmi_sb_t *)(dfmi->buf + wptr); wptr += nw_l;
    dfmi->sa_occ = (gin_dfmi_sao_t*)(dfmi->buf + wptr); wptr += nw_sa_occ;
    dfmi->sa     = (gin_uint40_t  *)(dfmi->buf + wptr); wptr += nw_sa;
    dfmi->unr    = (dfmi->buf + wptr);
    return dfmi;
}

int64_t gin_dfmi_bwt_length(gin_dfmi_t *dfmi) {
    return (int64_t)dfmi->header.values.no_chars-1;
}

void gin_dfmi_populate_alphabet(gin_dfmi_t *dfmi, int64_t** alphabet, int64_t *alphabet_size) {
    *alphabet_size = DFMI_NO_ALL_CHARS;
    *alphabet = (int64_t*) calloc(*alphabet_size, sizeof(int64_t));
    for(uint64_t i = 0; i < DFMI_NO_ALL_CHARS; i++) {
        (*alphabet)[i] = (int64_t)gin_dfmi_dec_lex[i];
    }
}

int64_t gin_dfmi_char_sa_base(gin_dfmi_t *dfmi, char c) {
    return (int64_t)dfmi->f[gin_dfmi_enc[c]];
}

int64_t gin_dfmi_size_in_bytes(gin_dfmi_t *dfmi) {
    return (int64_t)dfmi->header.values.fmi_size_in_bytes;
}

void gin_dfmi_decode(gin_dfmi_t *dfmi, char **string, uint64_t *len) {
    fprintf(stderr, "[gin_dna_fmi.c] Error: Decode for DFMI is not implemented. Quitting.\n");
    exit(EXIT_FAILURE);
}

int gin_dfmi_comp(gin_dfmi_t *dfmi1, gin_dfmi_t *dfmi2) {
    if (dfmi1->header.values.fmi_size_in_bytes == dfmi2->header.values.fmi_size_in_bytes) {
        return memcmp(dfmi1->buf, dfmi2->buf, dfmi1->header.values.fmi_size_in_bytes);
    }
    return -1;
}

uint64_t gin_dfmi_hash(gin_dfmi_t *dfmi) {
    const uint64_t prime = 1099511628211LLU;
    uint64_t hash = 14695981039346656037LLU;
    uint64_t *seq64 = dfmi->buf;
    for (uint64_t i = 0; i < dfmi->header.values.fmi_size_in_bytes / sizeof(uint64_t); ++i) {
        hash ^= seq64[i];
        hash *= prime;
    }
    return hash;
}

void gin_dfmi_free(gin_dfmi_t *dfmi) {
    free(dfmi->buf);
    free(dfmi);
}

void* gin_dfmi_copy(gin_dfmi_t *dfmi) {
    gin_dfmi_t *copy = calloc(1, sizeof(gin_dfmi_t));
    if(!copy) return NULL;
    uint8_t *buf = malloc(sizeof(uint8_t) * dfmi->header.values.fmi_size_in_bytes);
    // copy the header
    copy->header = dfmi->header;
    // copy the buffer
    memcpy(buf, dfmi->buf, dfmi->header.values.fmi_size_in_bytes);
    copy->buf = (uint64_t*)buf;
    // set pointers
    // Calculate sizes
    uint64_t N = dfmi->header.values.no_chars;
    uint64_t n_sb        = DFMI_CEIL(N,192);
    uint64_t n_sa_occ_sb = DFMI_CEIL(N,384);
    uint64_t n_sa_val =    DFMI_CEIL(N, dfmi->header.values.no_chars);
    // Number of words needed to store things
    uint64_t nw_header = sizeof(gin_dfmi_header_t) / sizeof(uint64_t);
    uint64_t nw_f = 8;
    uint64_t nw_l = n_sb * sizeof(gin_dfmi_sb_t) / sizeof(uint64_t);
    uint64_t nw_sa_occ = n_sa_occ_sb * sizeof(gin_dfmi_sao_t) / sizeof(uint64_t);
    uint64_t nw_sa = DFMI_CEIL(n_sa_val * sizeof(gin_uint40_t), sizeof(uint64_t));
    uint64_t wptr = nw_header;
    copy->f      = (copy->buf + wptr); wptr += nw_f;
    copy->l      = (gin_dfmi_sb_t *)(copy->buf + wptr); wptr += nw_l;
    copy->sa_occ = (gin_dfmi_sao_t*)(copy->buf + wptr); wptr += nw_sa_occ;
    copy->sa     = (gin_uint40_t  *)(copy->buf + wptr); wptr += nw_sa;
    copy->unr    = (copy->buf + wptr);
    return copy;
}

uint64_t gin_dfmi_count(gin_dfmi_t *dfmi, char *str) {
    uint64_t lo = 0;
    uint64_t hi = dfmi->header.values.no_chars+1;
    int64_t pos = (int64_t)strlen(str)-1;
    while(pos > -1) {
        char c = str[pos];
        uint8_t enc = gin_dfmi_enc[c];
        uint64_t rank_lo_m_1 = lo ? gin_dfmi_rank(dfmi,lo-1,c) : 0ull;
        uint64_t rank_hi_m_1 = hi ? gin_dfmi_rank(dfmi,hi-1,c) : 0ull;
        uint64_t base = dfmi->f[enc];
        lo = (base + rank_lo_m_1);
        hi = (base + rank_hi_m_1);
        if(hi <= lo) {
            lo = hi;
            break;
        }
        --pos;
    }
    return hi-lo;
}

void gin_dfmi_locate(gin_dfmi_t *dfmi, char *str, uint64_t **locs, uint64_t *nlocs) {
    uint64_t lo = 0;
    uint64_t hi = dfmi->header.values.no_chars+1;
    int64_t pos = (int64_t)strlen(str)-1;
    while(pos > -1) {
        char c = str[pos];
        uint8_t enc = gin_dfmi_enc[c];
        uint64_t rank_lo_m_1 = lo ? gin_dfmi_rank(dfmi,lo-1,c) : 0ull;
        uint64_t rank_hi_m_1 = hi ? gin_dfmi_rank(dfmi,hi-1,c) : 0ull;
        uint64_t base = dfmi->f[enc];
        lo = (base + rank_lo_m_1);
        hi = (base + rank_hi_m_1);
        if(hi <= lo) {
            lo = hi;
            break;
        }
        --pos;
    }
    if(hi > lo) {
        *locs = calloc(hi-lo, sizeof(uint64_t));
        *nlocs = hi - lo;
        gin_dfmi_sa(dfmi, *locs, lo, hi-1);
    } else {
        *locs = NULL;
        *nlocs = 0;
    }
}


void gin_dfmi_double_rank(gin_dfmi_t *dfmi, uint64_t pos, char c1, char c2, int64_t *r1, int64_t *r2) {
    // Encode chars
    uint8_t enc1 = gin_dfmi_enc[c1];
    uint8_t enc2 = gin_dfmi_enc[c2];
    // Then, go to the superblock and retrieve
    // a) superblock counts b) block counts c) the wavelets
    uint64_t s  =  pos / 192;
    uint64_t b  = (pos % 192) / 64;
    uint64_t m  =  pos % 64;
    // Character 1
    uint64_t sbc1 = DFMI_UINT40_GET(dfmi->l[s].sb.sbc[enc1]);
    uint64_t bc1  = (uint64_t)dfmi->l[s].sb.blocks[b].bc[enc1];
    uint64_t wav1 = gin_dfmi_wavelet_enc(enc1, dfmi->l[s].sb.blocks[b].bv);
    // Character 2
    uint64_t sbc2 = DFMI_UINT40_GET(dfmi->l[s].sb.sbc[enc2]);
    uint64_t bc2  = (uint64_t)dfmi->l[s].sb.blocks[b].bc[enc2];
    uint64_t wav2 = gin_dfmi_wavelet_enc(enc2, dfmi->l[s].sb.blocks[b].bv);
    // Popcounts
    uint64_t pc1 = DFMI_POPCOUNT(wav1 & DFMI_POPCOUNT_MASK(m));
    uint64_t pc2 = DFMI_POPCOUNT(wav2 & DFMI_POPCOUNT_MASK(m));
    // Mask the wavelet up until pos, and return the sum of cached counts
    *r1 = (int64_t)(sbc1 + bc1 + pc1);
    *r2 = (int64_t)(sbc2 + bc2 + pc2);
}
