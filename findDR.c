#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"
#include <time.h>

/******************************
 *  Direct Repeat Finder   ****
 *****************************/
int findDR(int mindir, int maxdir, int dspacer, int total_bases) {

	time_t rawtime;
	time(&rawtime);

	register int j, i, k, sp, end;
	int strti = 0;
	int ndx = 0;
	int size = 0;
	int lasti;
	int sizeMin = 0;
	int spMin;
	int spMax;
	int totlen;
	i = j = k = sp = size = end = 0;

	/*******************************************
	 * Start looking for direct repeats ********
	 *******************************************
	 /*/
	lasti = total_bases - (mindir * 2); //last i value that needs to be examined

	for (strti = 0; strti <= lasti; strti++) {
		while (dna[strti] == 'n') {//skip n's
			strti++;
		}
		if (strti >= lasti) {
			break;
		}

		for (size = maxdir; size >= mindir; size--) {
			if (((size * 2) + dspacer) <= (end - strti)) {// sizes are not big enough to escape prev.
				sp = dspacer;
				size = sizeMin;
				continue;
				//break;
			}
			//only examine spacers that give large enough repeats to escape previous
			spMin = max(0,((end-strti)-(size*2))+2);
			//watch for end of sequence
			spMax = min(dspacer,lasti-strti);
			for (sp = spMin; sp <= spMax; sp++) {
				j = strti + size + sp;
				i = strti;
				k = 0;
				while (dna[i] == dna[j] && k < size && dna[i] != 'n' && j
						< total_bases) { //using while so we can not compare all if not required
					k++;
					j++;
					i++;
				}
				if (k == size) {//DR found!
					totlen = k;
					if (sp == 0) {
						while (dna[i] == dna[j]) {//expand Big Tandem Repeat (BTR)
							totlen++;
							j++;
							i++;
						}
					}
					//Set new DR
					drep[ndx].start = strti + 1;
					drep[ndx].len = size;
					drep[ndx].loop = sp;
					drep[ndx].num = totlen / drep[ndx].len; //repeats
					drep[ndx].end = j;
					drep[ndx].sub = (totlen % drep[ndx].len); //remainder
					drep[ndx].strand = 0;

					ndx++;
					end = j - 1;
					sp = dspacer;
					size = sizeMin;
				}
			}
		}
	}
	return (ndx);
}/* END of direct*/
