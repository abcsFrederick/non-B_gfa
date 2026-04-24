#include <stdio.h>
#include <stdlib.h>
#include "gfa.h"

static void delIRep(REP *irep, int nreps, int toRemove) {
    int i;
    for (i = toRemove; i < nreps; i++)
        irep[i] = irep[i + 1];
}

int findIR(const char *dna, const char *dna3, REP *irep,
           int mincrf, int cspacer, int cut, int shortSpacer,
           int total_bases) {

    register int i, j, k, sp;
    int strti = 0;
    int cBack = 0;
    int maxcBack = 10;
    int ndx = 0;
    int tmpStart = 0;
    int tmpStop  = 0;
    int maxSP    = 0;
    i = j = k = sp = 0;
    BOOLEAN rightShifted = FALSE;
    BOOLEAN leftShifted  = FALSE;

    for (strti = mincrf; strti <= (total_bases - mincrf); strti++) {
        maxSP = min(cspacer, (total_bases - (strti + mincrf)));
        for (sp = 0; sp <= maxSP; sp++) {
            i = strti;
            k = 0;
            j = strti + sp + 1;
            while ((dna[i] == dna3[j]) && (j < total_bases) && (i >= 0)
                   && (dna[j] != 'n')) {
                k++;
                j++;
                i--;
            }
            if (k >= mincrf) {
                if ((k <= cut) && (sp > shortSpacer))
                    continue;

                tmpStart = ((strti - k) + 2);
                tmpStop  = (strti + k + sp + 1);

                if (ndx == 0) {
                    rightShifted = FALSE;
                    leftShifted  = FALSE;
                    irep[ndx].start  = tmpStart;
                    irep[ndx].sub    = tmpStop;
                    irep[ndx].len    = k;
                    irep[ndx].loop   = sp;
                    irep[ndx].num    = 1;
                    irep[ndx].end    = tmpStop;
                    irep[ndx].strand = 0;
                    ndx++;
                } else {
                    while (ndx > 0
                           && (irep[ndx-1].end <= tmpStop)
                           && (irep[ndx-1].start >= tmpStart)
                           && irep[ndx-1].len < k) {
                        ndx--;
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    }
                    while (ndx > 0
                           && (irep[ndx-1].end >= tmpStop)
                           && (irep[ndx-1].start <= tmpStart)
                           && irep[ndx-1].len < k) {
                        ndx--;
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    }
                    if (ndx > 0
                        && (irep[ndx-1].end <= tmpStop)
                        && (irep[ndx-1].start >= tmpStart)
                        && irep[ndx-1].len > k) {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    } else if (ndx > 0
                               && (irep[ndx-1].end >= tmpStop)
                               && (irep[ndx-1].start <= tmpStart)
                               && irep[ndx-1].len > k) {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    } else if (ndx > 0
                               && (tmpStop == irep[ndx-1].end)
                               && (k == irep[ndx-1].len)
                               && (!rightShifted)) {
                        leftShifted  = TRUE;
                        rightShifted = FALSE;
                        ndx--;
                        irep[ndx].num++;
                        irep[ndx].sub = tmpStart;
                        ndx++;
                    } else if (ndx > 0
                               && (tmpStart == irep[ndx-1].start)
                               && (k == irep[ndx-1].len)
                               && (!leftShifted)) {
                        rightShifted = TRUE;
                        leftShifted  = FALSE;
                        if (ndx > 1
                            && (irep[ndx-2].end <= tmpStop)
                            && (irep[ndx-2].start >= tmpStart)
                            && irep[ndx-2].len < k) {
                            irep[ndx-2] = irep[ndx-1];
                            ndx--;
                        }
                        ndx--;
                        irep[ndx].num++;
                        irep[ndx].end = tmpStop;
                        ndx++;
                    } else {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                        irep[ndx].start  = tmpStart;
                        irep[ndx].sub    = tmpStop;
                        irep[ndx].len    = k;
                        irep[ndx].loop   = sp;
                        irep[ndx].num    = 1;
                        irep[ndx].end    = tmpStop;
                        irep[ndx].strand = 0;
                        ndx++;

                        for (cBack = 1; cBack <= maxcBack; cBack++) {
                            while (((cBack+1) <= ndx)
                                   && (((irep[ndx-(1+cBack)].end >= irep[ndx-1].end)
                                        && (irep[ndx-(1+cBack)].start <= irep[ndx-1].start))
                                       || ((irep[ndx-(1+cBack)].end <= irep[ndx-1].end)
                                           && (irep[ndx-(1+cBack)].start >= irep[ndx-1].start)))) {
                                if (irep[ndx-(1+cBack)].len == irep[ndx-1].len) {
                                    if (irep[ndx-(1+cBack)].loop > irep[ndx-1].loop)
                                        delIRep(irep, (ndx-1), (ndx-(1+cBack)));
                                } else if (irep[ndx-(1+cBack)].len < irep[ndx-1].len) {
                                    delIRep(irep, (ndx-1), (ndx-(1+cBack)));
                                }
                                rightShifted = FALSE;
                                leftShifted  = FALSE;
                                --ndx;
                                cBack = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return ndx;
}
