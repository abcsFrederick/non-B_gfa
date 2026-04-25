#include <stdio.h>
#include "gfa.h"

/* Compute the reverse-complement strand (dna → dna2). */
void rcdna(const char *dna, char *dna2, int ndna) {
    int i, k;
    for (i = 0; i < ndna; i++) {
        k = ndna - i - 1;
        switch (dna[i]) {
            case 'a': dna2[k] = 't'; break;
            case 'c': dna2[k] = 'g'; break;
            case 'g': dna2[k] = 'c'; break;
            case 't': dna2[k] = 'a'; break;
            default:  dna2[k] = 'n'; break;
        }
    }
    dna2[ndna] = '\0';
}
