#include <stdio.h>
#include <stdlib.h>
#include "gfa.h"

static int pupy(const char *dna, int pos) {
    int isPPY = 0;
    if (dna[pos] == 'a') {
        if (dna[pos+1] == 'c') isPPY = 3;
    } else if (dna[pos] == 't') {
        if (dna[pos+1] == 'g') isPPY = 3;
    } else if (dna[pos] == 'c') {
        if      (dna[pos+1] == 'g') isPPY = 25;
        else if (dna[pos+1] == 'a') isPPY = 3;
    } else if (dna[pos] == 'g') {
        if      (dna[pos+1] == 'c') isPPY = 25;
        else if (dna[pos+1] == 't') isPPY = 3;
    }
    return isPPY;
}

int findZDNA(const char *dna, REP *zrep, int minZ, int total_bases) {
    register int i;
    int ndx    = 0;
    int tmpPPY = 0;
    int npy    = 1;
    int kvsum  = 0;

    i = 0;
    while (i < (total_bases - minZ)) {
        tmpPPY = pupy(dna, i);
        if (tmpPPY > 0) {
            npy++;
            kvsum += tmpPPY;
        } else {
            if (npy >= minZ) {
                zrep[ndx].start  = i - npy + 2;
                zrep[ndx].len    = npy;
                zrep[ndx].loop   = kvsum / 2;
                zrep[ndx].num    = 0;
                zrep[ndx].end    = i + 1;
                zrep[ndx].sub    = 0;
                zrep[ndx].strand = 0;
                ndx++;
            }
            npy   = 1;
            kvsum = 0;
        }
        i++;
    }
    return ndx;
}
