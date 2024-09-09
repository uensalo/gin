#include "gin_dna_fmi.h"
#include "divsufsort64.h"

gin_dfmi_t* gin_dfmi_build(const char* str, uint64_t size, uint64_t isa_rate) {
    gin_dfmi_t *dfmi = calloc(1, sizeof(gin_dfmi_t));
    if(!dfmi) return NULL;
    uint64_t nb = ((size-1) >> 8) + 1;
    uint64_t nsb = ((nb-1) >> 1) + 1;
    uint64_t nbv_words = ((size-1) >> 6) + 1;
    uint64_t n_isa = (size-1) / isa_rate + 1;
    uint64_t n_sa_occ_bv_blocks = ((size-1) >> 9) + 1;
    uint64_t n_sa_occ_bv_words = nbv_words;
    uint64_t n_total_words = 10 +                   // no_chars, isa_r, F
                             8 * nsb +              // headers, L
                             6 * nbv_words +        // bitvectors, L
                             n_sa_occ_bv_blocks +   // popcount for sa occupancy (1 word per block)
                             n_sa_occ_bv_words +    // sa occupancy bitvector
                             n_isa +                // sa entries
                             nbv_words;             // permutation
    uint64_t *buf = calloc(n_total_words, sizeof(uint64_t));

    // populate starts and write simple stuff to the bs
    dfmi->buf = buf;
    buf[0] = size;
    buf[1] = isa_rate;
    dfmi->no_chars = size;
    dfmi->isa_r = isa_rate;

    //pointers to places
    dfmi->f = &buf[2];
    dfmi->headers = (gin_dfmi_header_t*)(&buf[10]);
    for(uint64_t i = 0; i < 4; i++) {
        dfmi->bv_1[i] = &buf[10 + (8 * nsb) + i * nbv_words];
    }
    for(uint64_t i = 0; i < 2; i++) {
        dfmi->bv_2[i] = &buf[10 + (8 * nsb) + (4 + i) * nbv_words];
    }
    dfmi->bv_sa_h   = &buf[10 + (8 * nsb) + (6 * nbv_words)];
    dfmi->bv_sa_occ = &buf[10 + (8 * nsb) + (6 * nbv_words) + n_sa_occ_bv_blocks];
    dfmi->sa   = &buf[10 + (8 * nsb) + (6 * nbv_words) + n_sa_occ_bv_blocks + n_sa_occ_bv_words];
    dfmi->bv_p = &buf[10 + (8 * nsb) + (7 * nbv_words) + n_sa_occ_bv_blocks + n_sa_occ_bv_words];
    // Compute the SA
    int64_t *sa = calloc(size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)str, (saidx64_t*)sa, (saidx64_t)size+1);

    char *bwt = (char*)malloc((size+1)*sizeof(char));
    // Compute the BWT
    uint64_t k = 0;
    for(uint64_t i = 0; i < size+1; i++) {
        bwt[i] = sa[i] ? str[sa[i] - 1] : (char)'\0';
        if(!(sa[i] % isa_rate)) {
            DFMI_SET_BIT(dfmi->bv_sa_occ, i);
            dfmi->sa[k++] = sa[i];
        }
    }
    // SA is no longer needed, free here
    free(sa);
    // Populate bitvectors
    for(uint64_t i = 0; i < size; i++) {

    }
}


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
