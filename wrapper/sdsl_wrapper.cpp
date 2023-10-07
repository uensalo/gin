#include "sdsl_wrapper.h"
#include <sdsl/suffix_arrays.hpp>


struct csa_wt_wrapper {
    sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>> fmi;
};

void* csa_wt_build(const char* str) {
    csa_wt_wrapper *wrapper = new csa_wt_wrapper();
    std::string input(str);
    sdsl::construct(wrapper->fmi, input, 1);
    return wrapper;
}

uint64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c) {
    csa_wt_wrapper *wrapper = static_cast<csa_wt_wrapper*>(obj_handle);
    // Call some function on the csa_wt object. Example:
    return wrapper->fmi.wavelet_tree.rank(pos, c); // Replace with actual function
}

void csa_wt_to_buffer(void* obj_handle, uint8_t *data, uint64_t size) {

}

void csa_wt_from_buffer(void* obj_handle, uint8_t *data, uint64_t size) {

}