#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gfa.h"

void getGislands(const char *dna,
                 G_Island *gisle,  int *nGisls,
                 G_Island *rcgisle, int *nCisls,
                 int minGQ, int total_bases) {
    register int i;
    int ngs = 0;
    int ncs = 0;

    *nGisls = 0;
    *nCisls = 0;
    i = 0;

    while (i <= total_bases) {
        if (dna[i] == 'g') {
            ngs++;
        } else {
            if (ngs >= minGQ) {
                gisle[*nGisls].strt = i - ngs + 1;
                gisle[*nGisls].len  = ngs;
                (*nGisls)++;
            }
            ngs = 0;
        }
        if (dna[i] == 'c') {
            ncs++;
        } else {
            if (ncs >= minGQ) {
                rcgisle[*nCisls].strt = i - ncs + 1;
                rcgisle[*nCisls].len  = ncs;
                (*nCisls)++;
            }
            ncs = 0;
        }
        i++;
    }
}

int findGQ(G_Island *gisle,  int nGisls,
           G_Island *rcgisle, int nCisls,
           REP *grep, int minGQ, int maxGQspacer) {
    int ndx = 0;
    int nIls = 0;
    G_Island *islands;
    register int i, i2;
    int npos, j, k, m;
    i = j = k = m = i2 = 0;
    int nposMax, maxGQ;
    nposMax = maxGQ = 0;
    npos = 0;
    int conIls;
    int strand;

    for (strand = 0; strand < 2; strand++) {
        if (strand == 0) { nIls = nGisls;  islands = &gisle[0]; }
        else             { nIls = nCisls;  islands = &rcgisle[0]; }

        for (i = 0; i < nIls; i++) {
            conIls = 1;
            npos = (int)(floor((islands[i].len + 1) / (minGQ + 1)));
            i2 = i + 1;
            while (((islands[i2].strt - (islands[i2-1].strt + islands[i2-1].len))
                    <= maxGQspacer) && (i2 < nIls)) {
                conIls++;
                npos += (int)(floor((islands[i2].len + 1) / (minGQ + 1)));
                i2++;
            }
            if (npos >= 4) {
                maxGQ = minGQ;
                for (j = i; j < i2; j++) {
                    for (k = islands[j].len; k > maxGQ; k--) {
                        nposMax = (int)(floor((islands[j].len + 1) / (k + 1)));
                        for (m = j + 1; m < i2; m++) {
                            nposMax += (int)(floor((islands[m].len + 1) / (k + 1)));
                            if (nposMax >= 4) {
                                maxGQ = k;
                                break;
                            }
                            if ((int)(floor((islands[m].len + 1) / (k + 1))) == 0) {
                                if (islands[m+1].strt > (islands[m-1].strt
                                        + islands[m-1].len + maxGQspacer))
                                    break;
                            }
                        }
                    }
                }
                grep[ndx].start  = islands[i].strt;
                grep[ndx].num    = npos;
                grep[ndx].sub    = conIls;
                grep[ndx].len    = maxGQ;
                grep[ndx].end    = (islands[i2-1].strt + islands[i2-1].len) - 1;
                grep[ndx].strand = strand;
                ndx++;
            }
            i = i + conIls - 1;
        }
    }
    return ndx;
}
