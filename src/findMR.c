#include <stdio.h>
#include <stdlib.h>
#include "gfa.h"

static void delMRep(REP *mrep, int nreps, int toRemove) {
    int i;
    for (i = toRemove; i < nreps; i++)
        mrep[i] = mrep[i + 1];
}

int findMR(const char *dna, REP *mrep,
           int minmir, int mspacer, int total_bases) {

    register int i, j, k, sp;
    int strti = 0;
    int cBack = 0;
    int maxcBack = 5;
    int ndx = 0;
    int tmpStart = 0;
    int tmpStop  = 0;
    int maxSP    = 0;
    i = j = k = sp = 0;
    BOOLEAN rightShifted = FALSE;
    BOOLEAN leftShifted  = FALSE;

    for (strti = minmir; strti <= (total_bases - minmir); strti++) {
        maxSP = min(mspacer, (total_bases - (strti + minmir)));
        for (sp = 0; sp <= maxSP; sp++) {
            i = strti;
            k = 0;
            j = strti + sp + 1;
            while ((dna[i] == dna[j]) && (j < total_bases) && (i >= 0)
                   && (dna[j] != 'n')) {
                k++;
                j++;
                i--;
            }
            if (k >= minmir) {
                tmpStart = ((strti - k) + 2);
                tmpStop  = (strti + k + sp + 1);

                if (ndx == 0) {
                    rightShifted = FALSE;
                    leftShifted  = FALSE;
                    mrep[ndx].start  = tmpStart;
                    mrep[ndx].sub    = tmpStop;
                    mrep[ndx].len    = k;
                    mrep[ndx].loop   = sp;
                    mrep[ndx].num    = 1;
                    mrep[ndx].end    = tmpStop;
                    mrep[ndx].strand = 0;
                    ndx++;
                } else {
                    while (ndx > 0
                           && (mrep[ndx-1].end <= tmpStop)
                           && (mrep[ndx-1].start >= tmpStart)
                           && mrep[ndx-1].len < k) {
                        ndx--;
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    }
                    while (ndx > 0
                           && (mrep[ndx-1].end >= tmpStop)
                           && (mrep[ndx-1].start <= tmpStart)
                           && mrep[ndx-1].len < k) {
                        ndx--;
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    }
                    if (ndx > 0
                        && (mrep[ndx-1].end <= tmpStop)
                        && (mrep[ndx-1].start >= tmpStart)
                        && mrep[ndx-1].len > k) {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    } else if (ndx > 0
                               && (mrep[ndx-1].end >= tmpStop)
                               && (mrep[ndx-1].start <= tmpStart)
                               && mrep[ndx-1].len > k) {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                    } else if (ndx > 0
                               && (tmpStop == mrep[ndx-1].end)
                               && (k == mrep[ndx-1].len)
                               && (!rightShifted)) {
                        leftShifted  = TRUE;
                        rightShifted = FALSE;
                        ndx--;
                        mrep[ndx].num++;
                        mrep[ndx].sub = tmpStart;
                        ndx++;
                    } else if (ndx > 0
                               && (tmpStart == mrep[ndx-1].start)
                               && (k == mrep[ndx-1].len)
                               && (!leftShifted)) {
                        rightShifted = TRUE;
                        leftShifted  = FALSE;
                        if (ndx > 1
                            && (mrep[ndx-2].end <= tmpStop)
                            && (mrep[ndx-2].start >= tmpStart)
                            && mrep[ndx-2].len < k) {
                            mrep[ndx-2] = mrep[ndx-1];
                            ndx--;
                        }
                        ndx--;
                        mrep[ndx].num++;
                        mrep[ndx].end = tmpStop;
                        ndx++;
                    } else {
                        rightShifted = FALSE;
                        leftShifted  = FALSE;
                        mrep[ndx].start  = tmpStart;
                        mrep[ndx].sub    = tmpStop;
                        mrep[ndx].len    = k;
                        mrep[ndx].loop   = sp;
                        mrep[ndx].num    = 1;
                        mrep[ndx].end    = tmpStop;
                        mrep[ndx].strand = 0;
                        ndx++;

                        for (cBack = 1; cBack <= maxcBack; cBack++) {
                            while (((cBack+1) <= ndx)
                                   && (((mrep[ndx-(1+cBack)].end >= mrep[ndx-1].end)
                                        && (mrep[ndx-(1+cBack)].start <= mrep[ndx-1].start))
                                       || ((mrep[ndx-(1+cBack)].end <= mrep[ndx-1].end)
                                           && (mrep[ndx-(1+cBack)].start >= mrep[ndx-1].start)))) {
                                if (mrep[ndx-(1+cBack)].len == mrep[ndx-1].len) {
                                    if (mrep[ndx-(1+cBack)].loop > mrep[ndx-1].loop)
                                        delMRep(mrep, (ndx-1), (ndx-(1+cBack)));
                                } else if (mrep[ndx-(1+cBack)].len < mrep[ndx-1].len) {
                                    delMRep(mrep, (ndx-1), (ndx-(1+cBack)));
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
