#include <stdio.h>
#include <stdlib.h>
#include "gfa.h"

int getAtracts(const char *dna, const char *dna2,
               potential_Bent_DNA *pAPRs,
               int minAT, int maxAT, int total_bases) {
    int n = 0, n_rc = 0;
    int nPATs = 0;

    int maxATlen = 0, maxTlen = 0;
    int Alen = 0, Tlen = 0, ATlen = 0, TAlen = 0, maxATend = 0, ATend = 0;

    int maxATlen_rc = 0, maxTlen_rc = 0;
    int Alen_rc = 0, Tlen_rc = 0, ATlen_rc = 0, TAlen_rc = 0, maxATend_rc = 0;

    register int i = 0;
    int nAs  = 0;
    int strt = 0;

    while (i < total_bases) {
        if ((dna[i] == 'a') || (dna[i] == 't')) {
            nAs++;
        } else {
            if ((nAs >= minAT) && (nAs <= maxAT)) {
                strt  = i - nAs + 1;
                ATend = strt + nAs;

                Alen = Tlen = ATlen = maxATlen = maxTlen = TAlen = maxATend = 0;
                Alen_rc = Tlen_rc = ATlen_rc = maxATlen_rc = maxTlen_rc
                        = TAlen_rc = maxATend_rc = 0;

                n_rc = total_bases - ATend;
                for (n = strt - 1; n < ATend - 1; n++) {
                    n_rc++;
                    if (dna[n] == 'a') {
                        Tlen = 0; TAlen = 0;
                        if (dna[n-1] == 't') { Alen = 0; ATlen = 0; }
                        else { Alen++; ATlen++; }
                    }
                    if (dna2[n_rc] == 'a') {
                        Tlen_rc = 0; TAlen_rc = 0;
                        if (dna2[n_rc-1] == 't') { Alen_rc = 0; ATlen_rc = 0; }
                        else { Alen_rc++; ATlen_rc++; }
                    }
                    if (dna[n] == 't') {
                        if (TAlen < Alen) { TAlen++; ATlen++; }
                        else { Tlen++; TAlen = 0; ATlen = 0; Alen = 0; }
                    }
                    if (dna2[n_rc] == 't') {
                        if (TAlen_rc < Alen_rc) { TAlen_rc++; ATlen_rc++; }
                        else { Tlen_rc++; TAlen_rc = 0; ATlen_rc = 0; Alen_rc = 0; }
                    }
                    if (maxATlen < ATlen)       { maxATlen = ATlen; maxATend = n; }
                    if (maxTlen  < Tlen)          maxTlen  = Tlen;
                    if (maxATlen_rc < ATlen_rc)  { maxATlen_rc = ATlen_rc; maxATend_rc = n_rc; }
                    if (maxTlen_rc  < Tlen_rc)    maxTlen_rc   = Tlen_rc;
                }

                if (((maxATlen - maxTlen) >= minAT)
                    || ((maxATlen_rc - maxTlen_rc) >= minAT)) {
                    pAPRs[nPATs].end  = strt + nAs;
                    pAPRs[nPATs].strt = strt;
                    if ((maxATlen - maxTlen) >= (maxATlen_rc - maxTlen_rc))
                        pAPRs[nPATs].a_center = ((double)maxATend
                                - (((double)maxATlen - 1) / 2)) + 1;
                    else
                        pAPRs[nPATs].a_center = total_bases
                                - (((double)maxATend_rc
                                        - (((double)maxATlen_rc - 1) / 2)));
                    nPATs++;
                }
            }
            nAs = 0;
        }
        i++;
    }
    return nPATs;
}

int findAPR(const char *dna, const char *dna2,
            REP *arep, potential_Bent_DNA *pAPRs,
            int minAPR, int maxAPR, int minATracts, int total_bases) {
    register int i = 0;
    int nProcessedATs;
    int tracts = 1;
    double distToNext = 0;
    int ndx = 0;

    nProcessedATs = getAtracts(dna, dna2, pAPRs, minAPR, maxAPR, total_bases);

    for (i = 0; i < nProcessedATs - (minATracts + 1); i++) {
        distToNext = pAPRs[i+1].a_center - pAPRs[i].a_center;
        if ((distToNext <= 11.1) && (distToNext >= 9.9)) {
            tracts++;
        } else {
            if (tracts >= minATracts) {
                arep[ndx].start  = pAPRs[(i - tracts) + 1].strt;
                arep[ndx].loop   = 0;
                arep[ndx].num    = tracts;
                arep[ndx].strand = 0;
                arep[ndx].len    = tracts;
                arep[ndx].end    = pAPRs[i].end - 1;
                ndx++;
            }
            tracts = 1;
        }
    }
    return ndx;
}
