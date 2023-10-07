#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <vector>

typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>> csa_type;

extern "C" {

void* csa_wt_build(const char* str) {
    csa_type* csa = new csa_type();
    sdsl::construct_im(*csa, std::string(str));
    return static_cast<void*>(csa);
}

uint64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c) {
    if(!obj_handle) return -1; // or some error indicator
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    return csa->bwt.rank(pos, c);
}

void csa_wt_sa(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end) {
    if(!obj_handle) return;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    for(uint64_t j = start; j <= end; ++j) {
        buf[j-start] = (*csa)[j];
    }
}

void csa_wt_to_buffer(void* obj_handle, uint8_t **data, uint64_t *size) {
    if(!obj_handle || !data || !size) return;

    csa_type* csa = static_cast<csa_type*>(obj_handle);
    std::stringstream ss;
    sdsl::serialize(*csa, ss);

    std::string serialized_data = ss.str();
    *size = serialized_data.size();

    uint64_t padded_size = (1 + (((*size) - 1) >> 3)) << 3;

    // Allocate memory for the buffer
    *data = (uint8_t*) calloc(padded_size, sizeof(uint8_t));
    if (!*data) {
        *size = 0;
        return; // Memory allocation failed
    }

    // Copy data to the buffer
    std::copy(serialized_data.begin(), serialized_data.end(), *data);
}

void* csa_wt_from_buffer(uint8_t *data, uint64_t size) {
    if(!data) return nullptr;

    std::string buf(reinterpret_cast<char*>(data), size);
    std::stringstream ss(buf);

    csa_type* csa = new csa_type();
    sdsl::load(*csa, ss);
    return static_cast<void*>(csa);
}

size_t csa_wt_bwt_length(void* obj_handle) {
    if (!obj_handle) return 0;

    csa_type* csa = static_cast<csa_type*>(obj_handle);
    return csa->size(); // This returns the length of the BWT (and the original string).
}

void csa_populate_alphabet(void* obj_handle, int64_t ** alphabet, int64_t *alphabet_size) {
    if (!obj_handle || !alphabet || !alphabet_size) return;

    csa_type* csa = static_cast<csa_type*>(obj_handle);

    std::vector<int64_t> temp_alphabet;
    for (int64_t ch = 0; ch < 256; ++ch) { // Assuming extended ASCII for simplicity
        if (csa->char2comp[ch] != 0xFF) { // 0xFF indicates the character is not present
            temp_alphabet.push_back(ch);
        }
    }

    *alphabet_size = static_cast<int64_t>(temp_alphabet.size());

    // Allocate memory for the alphabet using calloc
    *alphabet = (int64_t*) calloc(*alphabet_size, sizeof(int64_t));
    if (!*alphabet) {
        *alphabet_size = 0;
        return;  // Memory allocation failed
    }

    // Populate the alphabet
    for (int64_t i = 0; i < *alphabet_size; ++i) {
        (*alphabet)[i] = temp_alphabet[i];
    }
}

void csa_wt_free(void* obj_handle) {
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    delete csa;
}

bool csa_wt_comp(void* obj_handle1, void* obj_handle2) {
    if (!obj_handle1 || !obj_handle2) return false;

    csa_type* csa1 = static_cast<csa_type*>(obj_handle1);
    csa_type* csa2 = static_cast<csa_type*>(obj_handle2);

    return csa1 == csa2;
}

size_t csa_wt_hash(void* obj_handle) {
    if (!obj_handle) return 0;

    csa_type* csa = static_cast<csa_type*>(obj_handle);
    std::stringstream ss;
    sdsl::serialize(*csa, ss);

    std::string serialized_data = ss.str();
    return std::hash<std::string>{}(serialized_data);
}

void* csa_wt_copy(void* obj_handle) {
    if (!obj_handle) return nullptr;

    csa_type* csa_src = static_cast<csa_type*>(obj_handle);
    csa_type* csa_dest = new csa_type(*csa_src);

    return static_cast<void*>(csa_dest);
}

} // extern "C"
