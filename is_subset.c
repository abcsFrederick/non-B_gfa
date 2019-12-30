#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "gfa.h"

/*****************************************************
 * sort repeats and remove overlaps ******************
 * only remove overlaps if smaller is ****************
 * centered within larger ****************************
 *****************************************************
 Defalut Values
 int minCruciformRep = 35;  //limit
 int maxCruciformSpacer = 3;  //max_loop
 int minTriplexYRpercent = 10;  //limit
 int maxTriplexSpacer = 8;  //max_loop
 int maxSlippedSpacer = 0;  //max_loop
 int minKVscore = 35; //limit
 */

void is_subset(int nreps, char X, int max_loop, int limit) {
	register int i;
	int j =0;
	REP *rep;
	int nY = 0;// Pyrimidine C/T count
	int nR = 0;// Purine A/G count
	int Ypercent = 0;

	//set pointer to correct repeat array
	if (X == 'M') {
		rep = &mrep[0];
	}
	else if (X == 'D') {
		rep = &drep[0];
	}
	else if (X == 'Z') {
		rep = &zrep[0];
	}
	else if (X == 'I') {
		rep = &irep[0];
	}
	else {
		fprintf(stderr, " FATAL Error in is_subset, exiting program\n");
		exit(20);
	}

	for (i = 0; i < nreps; i++) {//go through all reps
		rep[i].special = 0;
		if (X == 'M') {
			//count Y and R's
			nY = 0;
			nR = 0;
			for (j = rep[i].start; j <= rep[i].start+rep[i].len; j++) {
				switch (dna[j]) {
					case 'a':
					nR++;
						break;
					case 'c':
					nY++;
						break;
					case 'g':
					nR++;
						break;
					case 't':
					nY++;
						break;
				}
			}
			if (nY == 0){
				Ypercent = 0;
			}
			else if(nR == 0){
				Ypercent = 100;
			}
			else{
				Ypercent=((nY/nR)*100);
			}

			if ((rep[i].loop < max_loop )&& (Ypercent<=limit)) {
				rep[i].special = 1;
			}
		}
		else if (X == 'D') {
			if (rep[i].loop <= max_loop) {
				rep[i].special = 1;
			}
		}
		else if (X == 'Z') {
			if (rep[i].loop >= limit) {
				rep[i].special = 1;
			}
		}
		else if (X == 'I') {
			if ((rep[i].len >= limit) && (rep[i].loop <= max_loop)) {
				rep[i].special = 1;
			}
		}
	}
} /* END of is_subset*/

