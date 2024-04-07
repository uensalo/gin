#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <vector>
#include "../wrapper/sdsl_wrapper.h"

typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::hyb_vector<>>> csa_type;

extern "C" {

void* csa_wt_build(const char* str, uint64_t size) {
    csa_type* csa = new csa_type();
    sdsl::construct_im(*csa, str, 1);
    //printf("%ld\n", csa->size());
    return static_cast<void*>(csa);
}

int64_t csa_wt_rank(void* obj_handle, uint64_t pos, char c) {
    if(!obj_handle) return -1; // or some error indicator
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    return csa->wavelet_tree.rank(pos, c);
}

void csa_wt_sa(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end) {
    if(!obj_handle) return;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    for(uint64_t j = start; j <= end; ++j) {
        buf[j-start] = (*csa)[j];
    }
}

void csa_wt_bwt(void* obj_handle, uint64_t *buf, uint64_t start, uint64_t end) {
    if (!obj_handle || !buf) return;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    for (uint64_t i = start; i <= end; ++i) {
        buf[i - start] = csa->bwt[i];
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
    uint64_t padded_words = padded_size >> 3;

    // Allocate memory for the buffer
    *data = (uint8_t*) calloc(padded_words, sizeof(uint64_t));
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

int64_t csa_wt_bwt_length(void* obj_handle) {
    if (!obj_handle) return 0;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    //printf("%ld\n",csa->size());
    return csa->size(); // This returns the length of the BWT (and the original string).
}

void csa_wt_populate_alphabet(void* obj_handle, int64_t ** alphabet, int64_t *alphabet_size) {
    if (!obj_handle || !alphabet || !alphabet_size) return;

    csa_type* csa = static_cast<csa_type*>(obj_handle);

    // Create a temporary vector to store the unique characters present in the BWT
    std::vector<int64_t> temp_alphabet;

    // Iterate over the possible characters and check if they exist in the BWT using the char2comp array
    for (int64_t ch = 0; ch < csa->sigma; ++ch) {
        int64_t actual_char = csa->comp2char[ch];
        if (actual_char != 0xFF) {
            temp_alphabet.push_back(actual_char);
        }
    }

    *alphabet_size = static_cast<int64_t>(temp_alphabet.size());

    // Allocate memory for the alphabet
    *alphabet = (int64_t*) calloc(*alphabet_size, sizeof(int64_t));
    if (!*alphabet) {
        *alphabet_size = 0;
        return;  // Memory allocation failed
    }

    // Copy the characters from the temporary vector to the allocated array
    std::copy(temp_alphabet.begin(), temp_alphabet.end(), *alphabet);
}

int64_t csa_wt_char_sa_base(void* obj_handle, char c) {
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    return csa->C[csa->char2comp[c]];
}

int64_t csa_wt_size_in_bytes(void* obj_handle) {
    if (!obj_handle) return 0;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    return (int64_t)sdsl::size_in_bytes(*csa);
}

void csa_wt_free(void* obj_handle) {
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    delete csa;
}

int csa_wt_comp(void* obj_handle1, void* obj_handle2) {
    if (!obj_handle1 || !obj_handle2) return false;

    csa_type* csa1 = static_cast<csa_type*>(obj_handle1);
    csa_type* csa2 = static_cast<csa_type*>(obj_handle2);

    return csa1 == csa2;
}

uint64_t csa_wt_hash(void* obj_handle) {
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

void csa_wt_decode(void *obj_handle, char **str, uint64_t *len) {
    if (!obj_handle) return;
    csa_type* csa = static_cast<csa_type*>(obj_handle);
    uint64_t original_length = csa->size();
    char *decoded = (char*) calloc(original_length+1, sizeof(char));
    if (!decoded) {
        *str = NULL;
        *len = 0;
        return;
    }
    uint64_t pos = 0;
    for (int64_t i = original_length-2; i >= 0; --i) {
        char ch = (char)csa->bwt[pos];
        if(!ch) {
            decoded[original_length] = 0;
        } else {
            decoded[i] = ch;
        }
        pos = csa->C[csa->char2comp[ch]] + csa->wavelet_tree.rank(pos + 1, ch) - 1;
    }
    decoded[original_length] = '\0';
    *str = decoded;
    *len = original_length;
}

} // extern "C"
