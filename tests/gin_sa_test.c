#include <time.h>
#include "gin_gin.h"

void timestamp() {
    time_t ltime;
    ltime=time(NULL);
    char *str = asctime( localtime(&ltime) );
    str[strlen(str)-1]=0;
    printf("[%s]: ", str);
}

gin_string_t *generate_random_string(int length) {
    char charset[] = "AFPOKJGNQPALERGHAJnjdljkGHIsrilonolp;uNOfrmkjobjmsrOGKJM";
    int charsetlen = strlen(charset);
    // generate uniformly random labels of random length
    gin_string_t *str;
    gin_string_init(&str, length);
    for(int i = 0; i < length; i++) {
        str->seq[i]= charset[rand() % charsetlen];
    }
    str->size = length;
    return str;
}

int test(int length, int rank_sample_rate, int isa_sample_rate, int sup_print) {
    if(!sup_print) {
        printf("================================================================================\n");
        printf("TEST: length=%d, R=%d, ISA=%d \n",length,rank_sample_rate,isa_sample_rate);
        printf("================================================================================\n");
        timestamp();printf("Generating random string\n");
    }
    gin_string_t *random_string = generate_random_string(length);
    if(!sup_print) {
        timestamp();printf("Generated. Constructing SA\n");
    }
    int64_t *sa = calloc(random_string->size+1, sizeof(int64_t));
    divsufsort64((const unsigned char*)random_string->seq, sa, random_string->size+1);
    if(!sup_print) {
        timestamp();printf("Done. Constructing FMI\n");
    }
    gin_fmi_t *fmi;
    gin_fmi_init_with_sa(&fmi, random_string, sa, rank_sample_rate, isa_sample_rate);
    int passed = 1;

    gin_fmi_qr_t q;
    q.lo = 0;
    q.hi = fmi->no_chars;
    gin_vector_t *safmi = gin_fmi_sa(fmi, &q);

    for(int_t i = q.lo; i < q.hi; i++) {
        int_t a = (int_t)safmi->data[i-q.lo];
        int_t b = (int_t)sa[i];
        if(a != b) {
            passed = 0;
            break;
        }
    }

    if(!sup_print) {
        char *str = passed ? "Pass" : "Fail";
        timestamp();printf("Test Result: %s\n", str);
    }

    gin_vector_free(safmi);
    gin_fmi_free(fmi);
    gin_string_free(random_string);
    return passed;
}

int main() {
    srand(60);
    int_t lengths[] = {1024,2048,4096,8192,16384,65536};
    int_t rank_rates[] = {1,4,16,64,256};
    int_t sa_rates[] = {1,2,4,8,16,32,64,128,256,512,1024};

    int result = 0;

    for(int_t i = 0; i < sizeof(lengths) / sizeof(int_t); i++) {
        for(int_t j = 0; j < sizeof(rank_rates) / sizeof(int_t); j++) {
            for(int_t k = 0; k < sizeof(sa_rates) / sizeof(int_t); k++) {
                result += test(lengths[i], rank_rates[j], sa_rates[k], false);
            }
        }
    }

    printf("Tests passed: %lld/%lld\n", result, sizeof(lengths) / sizeof(int_t) * sizeof(rank_rates) / sizeof(int_t) *  sizeof(sa_rates) / sizeof(int_t));
}