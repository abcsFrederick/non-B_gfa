#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa.h"

static int compar4(const void *a, const void *b) {
    int i = ((REP *)a)->start - ((REP *)b)->start;
    if (i == 0) i = ((REP *)a)->len - ((REP *)b)->len;
    return i;
}

static void removeSTR(REP *srep, int nreps, int toRemove) {
    int i;
    for (i = toRemove; i < nreps; i++)
        srep[i] = srep[i + 1];
}

static int nonBstr(const char *dna, int start, int len) {
    int code = 0;
    int j = 0, i = 0;
    BOOLEAN isEven     = FALSE;
    BOOLEAN isSymetric = TRUE;
    BOOLEAN isPUPY     = TRUE;
    BOOLEAN isComp     = TRUE;

    if (len % 2 == 0) isEven = TRUE;
    j = start + len - 2;

    if (len >= 2) {
        for (i = 0; i <= (len / 2) - 1; i++) {
            if (dna[start + i - 1] != dna[j]) isSymetric = FALSE;
            if (isEven) {
                if      (dna[start+i-1] == 'a' && dna[j] != 't') isComp = FALSE;
                else if (dna[start+i-1] == 't' && dna[j] != 'a') isComp = FALSE;
                else if (dna[start+i-1] == 'c' && dna[j] != 'g') isComp = FALSE;
                else if (dna[start+i-1] == 'g' && dna[j] != 'c') isComp = FALSE;
            } else {
                isComp = FALSE;
            }
            j--;
        }
        for (i = start; i < (len + start - 1); i++) {
            if ((dna[i] == 'a') || (dna[i] == 'g')) {
                if ((dna[i-1] == 'a') || (dna[i-1] == 'g')) isPUPY = FALSE;
            }
            if ((dna[i] == 't') || (dna[i] == 'c')) {
                if ((dna[i-1] == 't') || (dna[i-1] == 'c')) isPUPY = FALSE;
            }
        }
    } else {
        isComp     = FALSE;
        isSymetric = TRUE;
        isPUPY     = FALSE;
    }

    if (isEven)     code += 1;
    if (isPUPY)     code += 2;
    if (isSymetric) code += 4;
    if (isComp)     code += 8;
    return code;
}

static int filterSTRs(REP *srep, int nSTRs) {
    int i;
    qsort(srep, nSTRs, sizeof(*srep), compar4);
    for (i = 1; i < nSTRs; i++) {
        if (srep[i].end <= srep[i-1].end) {
            removeSTR(srep, nSTRs, i);
            nSTRs--;
            i--;
        }
    }
    return nSTRs;
}

int findSTR(const char *dna, REP *srep,
            int minSTR, int maxSTR, int minSTRlen, int minReps,
            int total_bases) {
    register int i, j;
    j = i = 0;
    int ndx       = 0;
    int rpsz      = 0;
    int reps      = 1;
    int remainder = 0;
    int rs = 0, re = 0;

    for (i = 0; i < (total_bases - minSTRlen); i++) {
        while (dna[i] == 'n' && i != (total_bases - 1)) i++;
        for (rpsz = minSTR; rpsz <= maxSTR; rpsz++) {
            reps = 1;
            j = i + rpsz;
            while (strncmp(&dna[i], &dna[j], rpsz) == 0) {
                reps++;
                j += rpsz;
                if (j + rpsz >= total_bases) {
                    fprintf(stderr, "out of bounds 1\n");
                    break;
                }
            }
            if (reps >= minReps) {
                remainder = 0;
                rs = i;
                re = j;
                while (dna[rs] == dna[re]) { remainder++; rs++; re++; }
                if (((reps * rpsz) + remainder) >= minSTRlen) {
                    if (ndx >= 1) {
                        if (srep[ndx-1].end < re) {
                            srep[ndx].start  = i + 1;
                            srep[ndx].end    = re;
                            srep[ndx].num    = reps;
                            srep[ndx].loop   = nonBstr(dna, i + 1, rpsz);
                            srep[ndx].len    = rpsz;
                            srep[ndx].sub    = remainder;
                            srep[ndx].strand = 0;
                            ndx++;
                            i = re - minSTRlen + 1;
                            break;
                        }
                    } else {
                        srep[ndx].start  = i + 1;
                        srep[ndx].end    = re;
                        srep[ndx].num    = reps;
                        srep[ndx].loop   = nonBstr(dna, i + 1, rpsz);
                        srep[ndx].len    = rpsz;
                        srep[ndx].sub    = remainder;
                        srep[ndx].strand = 0;
                        ndx++;
                        i = re - minSTRlen + 1;
                        break;
                    }
                }
            }
        }
    }
    return filterSTRs(srep, ndx);
}
