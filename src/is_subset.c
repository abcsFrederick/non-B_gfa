#include <stdio.h>
#include <stdlib.h>
#include "gfa.h"

void is_subset(const char *dna, REP *rep, int nreps,
               char X, int max_loop, int limit) {
    register int i;
    int j = 0;
    int nY = 0, nR = 0, Ypercent = 0;

    for (i = 0; i < nreps; i++) {
        rep[i].special = 0;

        if (X == 'M') {
            nY = 0; nR = 0;
            for (j = rep[i].start; j <= rep[i].start + rep[i].len; j++) {
                switch (dna[j]) {
                    case 'a': nR++; break;
                    case 'c': nY++; break;
                    case 'g': nR++; break;
                    case 't': nY++; break;
                }
            }
            if      (nY == 0) Ypercent = 0;
            else if (nR == 0) Ypercent = 100;
            else              Ypercent = (nY / nR) * 100;

            if ((rep[i].loop < max_loop) && (Ypercent <= limit))
                rep[i].special = 1;

        } else if (X == 'D') {
            if (rep[i].loop <= max_loop)
                rep[i].special = 1;

        } else if (X == 'Z') {
            if (rep[i].loop >= limit)
                rep[i].special = 1;

        } else if (X == 'I') {
            if ((rep[i].len >= limit) && (rep[i].loop <= max_loop))
                rep[i].special = 1;

        } else {
            fprintf(stderr, "FATAL Error in is_subset: unknown type '%c'\n", X);
        }
    }
}
