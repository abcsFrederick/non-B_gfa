#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include "gfa.h"


/*****************************************************
 * comparStartEnd routine for qsort ******************
 * sort on the start then end position ***************
 *****************************************************
 */
int comparStartEnd(const void *a, const void *b) {
	int i = 0;
	i = ((REP *) a)->start - ((REP *) b)->start;
	if (i == 0) i = ((REP *) b)->end - ((REP *) a)->end;
	return (i);
}

/*****************************************************
 * remove a repeat, re-index rest of array************
 *****************************************************
 */
void removeRep(REP *rep, int nreps, int toRemove) {
	int i;
	for (i = toRemove; i < nreps; i++) {
		rep[i] = rep[i + 1];
	}
}

/*****************************************************
 * sort repeats and remove overlaps ******************
 * only remove overlaps if smaller is ****************
 * centered within larger ****************************
 *****************************************************
 */
int process_repeatsCentered(int nreps, char X) {
	register int i;
	REP *rep;

	//set pointer to correct repeat array
	if (X == 'M') {
		rep = &mrep[0];
	}
	if (X == 'I') {
		rep = &irep[0];
	}
	if (X == 'D') {
		rep = &drep[0];
	}

	// sort repeats by start then end locations
	qsort(rep, nreps, sizeof(*rep), comparStartEnd);

	//  nested if's to remove repeats which are part of
	for (i = 1; i < nreps; i++) {//go through all reps, start with 1 so i-1 == 0
		if (rep[i].end <= rep[i - 1].end) {//i within prev

			if ((rep[i].start - rep[i - 1].start) == (rep[i].end
					- rep[i - 1].end)) {//i is centered in prev?
				if (rep[i].start == rep[i - 1].start) {//complete overlap
					if (rep[i].loop >= rep[i - 1].loop) {//which rep has larger loop?  remove it
						removeRep(rep, nreps, i);
						nreps--;
						i--;//must re-check i-1 against next
					}
					else {//2: complete overlap, larger loop i-1
						removeRep(rep, nreps, i - 1);
						nreps--;
						i--;//must re-check i-1 against next
					}
				}
				else {// 3: i is centered subset, but not complete overlap
					removeRep(rep, nreps, i);
					nreps--;
					i--;//must re-check i-1 against next
				}
			}
			else { //i is un-centered subset of prev
				//annotate as sub/master
				rep[i].sub = 1;
				rep[i - 1].sub = 0;
			}
		}
	}
	return (nreps);
} /* END of process_repeatsCentered*/

//process centered and un-centered repeats within other repeats
int process_repeatsIncluded(int nreps, char X) {

	register int i;
	REP *rep;

	//set pointer to correct repeat array
	if (X == 'M') {
		rep = &mrep[0];
	}
	if (X == 'I') {
		rep = &irep[0];
	}
	if (X == 'D') {
		rep = &drep[0];
	}

	//  repeats sorted by end then start locations
	qsort(rep, nreps, sizeof(*rep), comparStartEnd);

	 //  nested if's to remove repeats which are part of larger repeats                                  **
	 for (i = 1; i < nreps; i++) {//go through all reps
		if (rep[i].end <= rep[i - 1].end) {//i within prev
			if ((rep[i].start - rep[i - 1].start) == (rep[i].end
					- rep[i - 1].end)) {//i is centered in prev?
				if (rep[i].start == rep[i - 1].start) {//complete overlap
					if (rep[i].loop >= rep[i - 1].loop) {//which rep has larger loop?  remove it
						//1: complete overlap, larger loop i
						removeRep(rep, nreps, i);
						nreps--;
						i--;//must re-check i-1 against next
					}
					else {//2: complete overlap, larger loop i-1
						removeRep(rep, nreps, i - 1);
						nreps--;
						i--;//must re-check i-1 against next
					}
				}//end if complete overlap
				else {// 3: i is centered subset, but not complete overlap
					removeRep(rep, nreps, i);
					nreps--;
					i--;//must re-check i-1 against next
				}
			}//end if centered
			else { //i is un-centered subset of prev, if prev is
				if (rep[i - 1].loop <= rep[i].loop) {//if loop i-1 less or equal to loop i
					//4: un-centered subset, larger stem is i-1
					removeRep(rep, nreps, i);
					nreps--;
					i--;//must re-check i-1 against next
				}
			}//end un-centered
		}
	}
	return (nreps);
} /* END of process_repeatsIncluded*/

