#ifndef FMD_SDSL_WRAPPER_H
#define FMD_SDSL_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>

void* csa_wt_build(const char* str);
uint64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c);
void csa_wt_to_buffer(void* obj_handle, uint8_t *data, uint64_t size);
void csa_wt_from_buffer(void* obj_handle, uint8_t *data, uint64_t size);

#ifdef __cplusplus
}
#endif

#endif //FMD_SDSL_WRAPPER_H
