#include <stdio.h>
#include <ctype.h>
#include "gfa.h"

/* Compute the complement strand (dna → dna3). */
void cdna(const char *dna, char *dna3, int ndna) {
    int i;
    for (i = 0; i < ndna; i++) {
        switch (dna[i]) {
            case 'a': dna3[i] = 't'; break;
            case 'c': dna3[i] = 'g'; break;
            case 'g': dna3[i] = 'c'; break;
            case 't': dna3[i] = 'a'; break;
            default:  dna3[i] = 'n'; break;
        }
    }
    dna3[ndna] = '\0';
}
