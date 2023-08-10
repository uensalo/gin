#include "fmd_fmi.h"

int main() {
    char *cstr = "AFPOKJGNQPALERGHAJnjdljkGHIsrilonolp;uNOfrmkjobjmsrOGKJM";
    fmd_string_t *str;
    fmd_string_init_cstr(&str, cstr);

    int64_t *sa = calloc(str->size+1, sizeof(uint64_t));
    divsufsort64((sauchar_t*)str->seq, (saidx64_t*)sa, str->size+1);

    for(int_t i = 0; i < str->size+1; i++) {
        printf("%ld ", sa[i]);
    }
    printf("\n");

    fmd_fmi_t *fmi;
    fmd_fmi_init_with_sa(&fmi, str, sa, 1, 1);

    fmd_fmi_qr_t q;
    q.lo = 0;
    q.hi = fmi->no_chars;
    fmd_vector_t *pos = fmd_fmi_sa(fmi, &q);
    for(int_t i = 0; i < pos->size; i++) {
        printf("%ld ", pos->data[i]);
    }
}