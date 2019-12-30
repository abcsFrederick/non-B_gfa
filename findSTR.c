#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "gfa.h"

/*******************************************
 * Start looking for Short Tandem Repeats **
 *******************************************
 */

/*****************************************************
 * compar4 routine for qsort *************************
 * sort on the repeat size (loop) then start**********
 *****************************************************
 */

//rep then start
int compar4(const void *a, const void *b) {
	int i = 0;
	i = ((REP *) a)->start - ((REP *) b)->start;
	if (i == 0) i = ((REP *) a)->len - ((REP *) b)->len;
	return (i);
}

void removeSTR(int nreps, int toRemove) {
	int i = 0;
	for (i = toRemove; i < nreps; i++) {
		srep[i] = srep[i + 1];
	}
}

//examines str sequence, returns int code for which non-B it can form
int nonBstr(int start, int len) {
	int code = 0;
	int j = 0;
	int i = 0;
	BOOLEAN isEven = FALSE; //even
	BOOLEAN isSymetric = TRUE; //symetic; ata, aagg...
	BOOLEAN isPUPY = TRUE; //acgtgca, tatat...
	BOOLEAN isComp = TRUE; //at, gc, aatt, (even only)
	if (len % 2 == 0) {
		isEven = TRUE;
	}
	j = start + len - 2;
	if (len >= 2) {
		for (i = 0; i <= (len / 2) - 1; i++) {
			if (dna[start + i - 1] != dna[j]) {
				isSymetric = FALSE;
			}
			if (isEven) {
				if (dna[start + i - 1] == 'a') {
					if (dna[j] != 't') {
						isComp = FALSE;
					}
				}
				else if (dna[start + i - 1] == 't') {
					if (dna[j] != 'a') {
						isComp = FALSE;
					}
				}
				else if (dna[start + i - 1] == 'c') {
					if (dna[j] != 'g') {
						isComp = FALSE;
					}
				}
				else if (dna[start + i - 1] == 'g') {
					if (dna[j] != 'c') {
						isComp = FALSE;
					}
				}
			}
			else {//odd
				isComp = FALSE;
			}
			j--;
		}//end for
		for (i = start; i < (len + start - 1); i++) {//isPUPY?
			if ((dna[i] == 'a') || (dna[i] == 'g')) {//R
				if ((dna[i - 1] == 'a') || (dna[i - 1] == 'g')) {//RR
					isPUPY = FALSE;
				}
			}
			if ((dna[i] == 't') || (dna[i] == 'c')) {//Y
				if ((dna[i - 1] == 't') || (dna[i - 1] == 'c')) {//YY
					isPUPY = FALSE;
				}
			}
		}
	}
	else {//length = 1
		isComp = FALSE;
		isSymetric = TRUE;
		isPUPY = FALSE;
	}

	//convert binary to decimal, first digit is even,
	// second is pupy, third is symetric, fourth is complement
	//0000 = 0, 0001 = 1, 0010 = 2, 0011 = 3, 0100 = 4, 0101 = 5,
	//0110 = 6, 0111 = 7, 1000 = 8, 1001 = 9, 1010 = 10, 1011 = 11,
	//1100 = 12, 1101 = 13, 1110 = 14, 1111 = 15

	if (isEven) {
		code = code + 1;
	}
	if (isPUPY) {
		code = code + 2;
	}
	if (isSymetric) {
		code = code + 4;
	}
	if (isComp) {
		code = code + 8;
	}
	return (code);
}

int filterSTRs(int nSTRs) {
	//sort by rep size then start
	qsort(srep, nSTRs, sizeof(*srep), compar4);
	//remove all that end before prev, will be
	int i = 0;
	for (i = 1; i < nSTRs; i++) {
		if (srep[i].end <= srep[i - 1].end) {
			removeSTR(nSTRs, i);
			nSTRs--;
			i--;
		}
	}

	return (nSTRs);
}

int findSTR(int minSTR, int maxSTR, int minSTRlen, int minReps, int total_bases) {
	register int i, j;
	j = 0;
	i = 0;
	int ndx = 0;
	int rpsz = 0;
	int reps = 1;
	int remainder = 0;
	int rs = 0;//remainder search start
	int re = 0;//remainder search end
	for (i = 0; i < (total_bases - minSTRlen); i++) {//for each nucleotide
		while (dna[i] == 'n' && i != (total_bases - 1)) {//skip any n stretches
			i++;
		}
		for (rpsz = minSTR; rpsz <= maxSTR; rpsz++) {//for each rep size
			reps = 1;
			j = i + rpsz;
			while (strncmp(&dna[i], &dna[j], rpsz) == 0) {
				reps++;
				j = j + rpsz;
				if (j + rpsz >= total_bases) {
					fprintf(stderr, "out of bounds 1\n");
					break;
				}
			}
			if (reps >= minReps) {
				remainder = 0;
				rs = i;
				re = j;
				while (dna[rs] == dna[re]) {//partial rep check
					remainder++;
					rs++;
					re++;
				}
				if ((((reps * rpsz) + remainder) >= minSTRlen)) {
					if (ndx >= 1) {//don't compare to previous if there isn't one
						if (srep[ndx - 1].end < re) {//make sure we escape last one
							srep[ndx].start = i + 1;//+1 because dna[] starts at 0
							srep[ndx].end = re;
							srep[ndx].num = reps;
							srep[ndx].loop = nonBstr(i + 1, rpsz);
							srep[ndx].len = rpsz;
							srep[ndx].sub = remainder;
							srep[ndx].strand = 0;
							ndx++;
							i = re - minSTRlen + 1;//skip to end of rep
							break;//for repsz
						}
					}
					else {
						srep[ndx].start = i + 1;//+1 because dna[] starts at 0
						srep[ndx].end = re;
						srep[ndx].num = reps;
						srep[ndx].loop = nonBstr(i + 1, rpsz);
						srep[ndx].len = rpsz;
						srep[ndx].sub = remainder;
						srep[ndx].strand = 0;
						ndx++;
						i = re - minSTRlen + 1;//skip to end of rep
						break;//for repsz
					}
				}
			}
		}
	}
	return (ndx);
}
