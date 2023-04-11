#ifndef FMD_FMD_COMMON_H
#define FMD_FMD_COMMON_H
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>

#define MIN2(a,b) ((a) <= (b) ? (a) : (b))
#define MIN3(a,b,c) (MIN2((MIN2((a),(b))),(c)))

#define MAX2(a,b) ((a) >= (b) ? (a) : (b))
#define MAX3(a,b,c) (MAX2((MAX2((a),(b))),(c)))
#define MAX4(a,b,c,d) (MAX2((MAX3((a),(b),(c))),(d)))

typedef uint8_t  byte_t;
typedef int64_t  pos_t;
typedef uint64_t upos_t;

typedef int    (*fcomp)(void*, void*);
typedef upos_t (*fhash)(void*);
typedef void   (*ffree)(void*);
typedef void*  (*fcopy)(void*);
typedef void   (*ftrav_kv)(void *key, void *value, void *p);

typedef struct fmd_fstruct {
    fcomp comp_f;
    fhash hash_f;
    ffree free_f;
    fcopy copy_f;
} fmd_fstruct_t;

// basic functions
static int prm_comp_f(void *a, void *b) {
    return a < b ? -1 : (a == b ? 0 : 1);
}
static upos_t prm_hash_f(void *a) {
    return (pos_t)a;
}
static void prm_free(void *a) {
}
static void* prm_copy(void *val) {
    return val;
}
static fmd_fstruct_t prm_fstruct = {prm_comp_f, prm_hash_f, prm_free, prm_copy};

// swap utility
static inline void fmd_swap(void **a, void **b) {
    void *temp = *a;
    *a = *b;
    *b = temp;
}

static int64_t fmd_popcount64(uint64_t x) {
    x -= ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + (x >> 2 & 0x3333333333333333ull);
    return (int64_t)(((x + (x >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56);
}

static int64_t fmd_ceil_log2(uint64_t x) {
    static const unsigned long long t[6] = {
            0xFFFFFFFF00000000ull,
            0x00000000FFFF0000ull,
            0x000000000000FF00ull,
            0x00000000000000F0ull,
            0x000000000000000Cull,
            0x0000000000000002ull
    };
    uint64_t y = (((x & (x - 1)) == 0) ? 0 : 1);
    uint64_t j = 32;
    for (uint64_t i = 0; i < 6; i++) {
        uint64_t k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }
    return (int64_t)y;
}

static bool fmd_test_miller_rabin(pos_t n, pos_t d) {
    pos_t a = 2 + rand() % (n - 4);
    pos_t x = 1;
    pos_t temp = d;
    while (temp > 0) {
        if (temp % 2 == 1) {
            x = (x * a) % n;
        }
        a = (a * a) % n;
        temp /= 2;
    }
    if (x == 1 || x == n - 1) {
        return true;
    }
    while (d != n - 1) {
        x = (x * x) % n;
        d *= 2;

        if (x == 1) {
            return false;
        }
        if (x == n - 1) {
            return true;
        }
    }
    return false;
}

static bool fmd_is_prime(pos_t n) {
    const pos_t no_iter = 5; // error = (1/4)^no_iter
    if (n <= 1 || n == 4) {
        return false;
    }
    if (n <= 3) {
        return true;
    }
    pos_t d = n - 1;
    while (d % 2 == 0) {
        d /= 2;
    }
    for (int i = 0; i < no_iter; i++) {
        if (!fmd_test_miller_rabin(n, d)) {
            return false;
        }
    }
    return true;
}

static pos_t fmd_next_prime_after(pos_t n) {
    if (n <= 1) {
        return 2;
    }
    pos_t prime = n % 2 == 0 ? n + 1 : n + 2;
    while (!fmd_is_prime(prime)) {
        prime += 2;
    }
    return prime;
}


#endif //FMD_FMD_COMMON_H
