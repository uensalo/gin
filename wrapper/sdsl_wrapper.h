#ifndef FMD_SDSL_WRAPPER_H
#define FMD_SDSL_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>

void* csa_wt_build(const char* str);
int64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c);
void csa_wt_sa(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end);
void csa_wt_to_buffer(void* obj_handle, uint8_t **data, uint64_t *size);
void* csa_wt_from_buffer(uint8_t *data, uint64_t size);
int64_t csa_wt_bwt_length(void* obj_handle);
void csa_wt_populate_alphabet(void* obj_handle, int64_t** alphabet, int64_t *alphabet_size);

int csa_wt_comp(void* obj_handle1, void* obj_handle2);
size_t csa_wt_hash(void* obj_handle);
void csa_wt_free(void* obj_handle);
void* csa_wt_copy(void* obj_handle);

#ifdef __cplusplus
}
#endif

#endif //FMD_SDSL_WRAPPER_H
