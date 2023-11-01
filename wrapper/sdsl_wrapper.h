#ifndef GIN_SDSL_WRAPPER_H
#define GIN_SDSL_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>

void* csa_wt_build(const char* str, uint64_t size);
int64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c);
void csa_wt_sa(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end);
void csa_wt_bwt(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end);
void csa_wt_to_buffer(void* obj_handle, uint8_t **data, uint64_t *size);
void* csa_wt_from_buffer(uint8_t *data, uint64_t size);
int64_t csa_wt_bwt_length(void* obj_handle);
void csa_wt_populate_alphabet(void* obj_handle, int64_t** alphabet, int64_t *alphabet_size);
int64_t csa_wt_char_sa_base(void* obj_handle, char c);
int64_t csa_wt_size_in_bytes(void* obj_handle);

int csa_wt_comp(void* obj_handle1, void* obj_handle2);
uint64_t csa_wt_hash(void* obj_handle);
void csa_wt_free(void* obj_handle);
void* csa_wt_copy(void* obj_handle);

#ifdef __cplusplus
}
#endif

#endif //GIN_SDSL_WRAPPER_H
