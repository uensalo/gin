#include "gin_dna_fmi.h"
#include "divsufsort64.h"

gin_dfmi_t* gin_dfmi_build(const char* str, uint64_t size, uint64_t isa_rate) {
    gin_dfmi_t *dfmi = calloc(1, sizeof(gin_dfmi_t));
    if(!dfmi) return NULL;

    return dfmi;
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
