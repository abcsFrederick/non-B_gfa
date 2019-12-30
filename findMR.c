#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

/******************************
 *  findMR Repeat Finder   ****
 *
 *****************************/

void delMRep(int nreps, int toRemove) {//shifts stack down, but does not reset ndx
	int i = 0;
	for (i = toRemove; i < nreps; i++) {
		mrep[i] = mrep[i + 1];
	}
	//fprintf(stderr, " removed start %d stop %d  \n", mrep[toRemove].start,mrep[toRemove].end);

}

int findMR(int minmir, int mspacer, int total_bases) {

	register int i, j, k, sp;
	int strti = 0;
	int cBack = 0;//counter to look at previous elements for overlap etc.
	int maxcBack = 5;//
	int ndx = 0;
	int tmpStart = 0;
	int tmpStop = 0;
	int maxSP = 0;
	i = j = k = sp = 0;
	BOOLEAN rightShifted = FALSE;
	BOOLEAN leftShifted = FALSE;

	/*******************************************
	 * Start looking for mirrors ************
	 *******************************************
	 */
	for (strti = minmir; strti <= (total_bases - minmir); strti++) {
		maxSP = min(mspacer,(total_bases-(strti+minmir)));
		for (sp = 0; sp <= maxSP; sp++) {
			i = strti;
			k = 0;
			j = strti + sp + 1;
			while ((dna[i] == dna[j]) && (j < (total_bases)) && (i >= 0)
					&& (dna[j] != 'n')) {
				k++;
				j++;
				i--;
			}
			if (k >= minmir) {
				tmpStart = ((strti - k) + 2);// in ncbi coordinates (+1 to array cords)
				tmpStop = (strti + k + sp + 1);// in ncbi coordinates (+1 to array cords)
				if ((ndx == 0)) {//first one, can't compare current to prev if prev doesn't exist
					rightShifted = FALSE;
					leftShifted = FALSE;
					mrep[ndx].start = tmpStart;
					mrep[ndx].sub = tmpStop;//min loop boundary set to end by default
					mrep[ndx].len = k;
					mrep[ndx].loop = sp;
					mrep[ndx].num = 1;
					mrep[ndx].end = tmpStop;
					mrep[ndx].strand = 0;

					ndx++;
				}
				else {//Not first one

					//check for immediate inclusions
					//old within new, new larger looped
					while ((mrep[ndx - 1].end <= tmpStop)
							&& (mrep[ndx - 1].start >= tmpStart) && mrep[ndx
							- 1].len < k) {
						//old within new, new better
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, new better
					while ((mrep[ndx - 1].end >= tmpStop)
							&& (mrep[ndx - 1].start <= tmpStart) && mrep[ndx
							- 1].len < k) {
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;

					}
					//old within new, old better
					if ((mrep[ndx - 1].end <= tmpStop) && (mrep[ndx - 1].start
							>= tmpStart) && mrep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, old better
					if ((mrep[ndx - 1].end >= tmpStop) && (mrep[ndx - 1].start
							<= tmpStart) && mrep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					else if ((tmpStop == mrep[ndx - 1].end) && (k == mrep[ndx
							- 1].len) && (!rightShifted)) {
						leftShifted = TRUE;//to check for alternate shifting
						rightShifted = FALSE;
						ndx--;
						mrep[ndx].num = mrep[ndx].num + 1;//add one to Permutation count
						mrep[ndx].sub = tmpStart;//adjust minimum loop boundary
						//mrep[ndx].sub = -11;
						ndx++;
					}
					else if ((tmpStart == mrep[ndx - 1].start) && (k
							== mrep[ndx - 1].len) && (!leftShifted)) {
						rightShifted = TRUE;//to check for alternate shifting
						leftShifted = FALSE;

						//need check as rightShfited grows that it doesn't swallow previous
						if ((mrep[ndx - 2].end <= tmpStop)
								&& (mrep[ndx - 2].start >= tmpStart)
								&& mrep[ndx - 2].len < k) {
							//old within new, new better
							//replace old with current shifting
							mrep[ndx - 2] = mrep[ndx - 1];
							ndx--;//replace previous
						}

						ndx--;
						mrep[ndx].num = (mrep[ndx].num + 1);//add one to Permutation count
						mrep[ndx].end = tmpStop;//adjust end
						//mrep[ndx].sub = -12;
						ndx++;
					}
					else {//neither shifted, add new repeat
						rightShifted = FALSE;
						leftShifted = FALSE;
						mrep[ndx].start = tmpStart;
						mrep[ndx].sub = tmpStop;//min loop boundary set to end by default
						mrep[ndx].len = k;
						mrep[ndx].loop = sp;
						mrep[ndx].num = 1;
						mrep[ndx].end = tmpStop;
						mrep[ndx].strand = 0;

						ndx++;

						for (cBack = 1; cBack <= maxcBack; cBack++) {
							while ((((mrep[ndx - (1 + cBack)].end >= mrep[ndx
									- 1].end) && (mrep[ndx - (1 + cBack)].start
									<= mrep[ndx - 1].start)) || ((mrep[ndx - (1
									+ cBack)].end <= mrep[ndx - 1].end)
									&& (mrep[ndx - (1 + cBack)].start
											>= mrep[ndx - 1].start)))
									&& ((cBack + 1) <= ndx)) {
								//maximize stem length, then minimize loop length
								if ((mrep[ndx - (1 + cBack)].len == mrep[ndx
										- 1].len)) {//if stems are equal, keep shortest loop
									//if previous loop is larger, delete it
									if ((mrep[ndx - (1 + cBack)].loop
											> mrep[ndx - 1].loop)) {
										delMRep((ndx - 1), (ndx - (1 + cBack)));
									}
								}
								//otherwise keep longest stem

								else if ((mrep[ndx - (1 + cBack)].len
										< mrep[ndx - 1].len)) {//previous stem is shorter = delete it
									delMRep((ndx - 1), (ndx - (1 + cBack)));
								}
								//reset everything and remove end ndx
								rightShifted = FALSE;
								leftShifted = FALSE;
								--ndx;
								cBack = 1;
							}
						}//while
					}//for cBack
				}//else first
			}//if ndx > maxcBack + 2
		}//if maxcBack>0
	}//if k>minmir
	return (ndx);
}/* END of findMR*/
