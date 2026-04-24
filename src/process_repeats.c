#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa.h"

static int comparStartEnd(const void *a, const void *b) {
    int i = ((REP *)a)->start - ((REP *)b)->start;
    if (i == 0) i = ((REP *)b)->end - ((REP *)a)->end;
    return i;
}

static void removeRep(REP *rep, int nreps, int toRemove) {
    int i;
    for (i = toRemove; i < nreps; i++)
        rep[i] = rep[i + 1];
}

int process_repeatsCentered(REP *rep, int nreps) {
    register int i;

    qsort(rep, nreps, sizeof(*rep), comparStartEnd);

    for (i = 1; i < nreps; i++) {
        if (rep[i].end <= rep[i-1].end) {
            if ((rep[i].start - rep[i-1].start) == (rep[i].end - rep[i-1].end)) {
                if (rep[i].start == rep[i-1].start) {
                    if (rep[i].loop >= rep[i-1].loop) {
                        removeRep(rep, nreps, i);
                        nreps--; i--;
                    } else {
                        removeRep(rep, nreps, i-1);
                        nreps--; i--;
                    }
                } else {
                    removeRep(rep, nreps, i);
                    nreps--; i--;
                }
            } else {
                rep[i].sub   = 1;
                rep[i-1].sub = 0;
            }
        }
    }
    return nreps;
}

int process_repeatsIncluded(REP *rep, int nreps) {
    register int i;

    qsort(rep, nreps, sizeof(*rep), comparStartEnd);

    for (i = 1; i < nreps; i++) {
        if (rep[i].end <= rep[i-1].end) {
            if ((rep[i].start - rep[i-1].start) == (rep[i].end - rep[i-1].end)) {
                if (rep[i].start == rep[i-1].start) {
                    if (rep[i].loop >= rep[i-1].loop) {
                        removeRep(rep, nreps, i);
                        nreps--; i--;
                    } else {
                        removeRep(rep, nreps, i-1);
                        nreps--; i--;
                    }
                } else {
                    removeRep(rep, nreps, i);
                    nreps--; i--;
                }
            } else {
                if (rep[i-1].loop <= rep[i].loop) {
                    removeRep(rep, nreps, i);
                    nreps--; i--;
                }
            }
        }
    }
    return nreps;
}
