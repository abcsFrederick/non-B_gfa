#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gfa.h"

int findDR(const char *dna, REP *drep,
           int mindir, int maxdir, int dspacer, int total_bases) {

    register int j, i, k, sp, end;
    int strti  = 0;
    int ndx    = 0;
    int size   = 0;
    int lasti;
    int sizeMin = 0;
    int spMin, spMax;
    int totlen;
    i = j = k = sp = size = end = 0;

    lasti = total_bases - (mindir * 2);

    for (strti = 0; strti <= lasti; strti++) {
        while (dna[strti] == 'n') {
            strti++;
        }
        if (strti >= lasti)
            break;

        for (size = maxdir; size >= mindir; size--) {
            if (((size * 2) + dspacer) <= (end - strti)) {
                sp   = dspacer;
                size = sizeMin;
                continue;
            }
            spMin = max(0, ((end - strti) - (size * 2)) + 2);
            spMax = min(dspacer, lasti - strti);
            for (sp = spMin; sp <= spMax; sp++) {
                j = strti + size + sp;
                i = strti;
                k = 0;
                while (dna[i] == dna[j] && k < size
                       && dna[i] != 'n' && j < total_bases) {
                    k++;
                    j++;
                    i++;
                }
                if (k == size) {
                    totlen = k;
                    if (sp == 0) {
                        while (dna[i] == dna[j]) {
                            totlen++;
                            j++;
                            i++;
                        }
                    }
                    drep[ndx].start  = strti + 1;
                    drep[ndx].len    = size;
                    drep[ndx].loop   = sp;
                    drep[ndx].num    = totlen / drep[ndx].len;
                    drep[ndx].end    = j;
                    drep[ndx].sub    = totlen % drep[ndx].len;
                    drep[ndx].strand = 0;
                    ndx++;
                    end  = j - 1;
                    sp   = dspacer;
                    size = sizeMin;
                }
            }
        }
    }
    return ndx;
}
