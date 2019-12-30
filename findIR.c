#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

/*******************************************************************
 *  findIR eXplorer (crux)                                      *
 *  Program to locate possible cruciforms in nucleic acid sequence *
 *******************************************************************/

void delIRep(int nreps, int toRemove) {//shifts stack down, but does not reset ndx
	int i = 0;
	for (i = toRemove; i < nreps; i++) {
		irep[i] = irep[i + 1];
	}
	//fprintf(stderr, " removed start %d stop %d  \n", mrep[toRemove].start,mrep[toRemove].end);

}

int findIR(int mincrf, int cspacer, int cut, int shortSpacer, int total_bases) {

	register int i, j, k, sp;
	int strti = 0;
	int cBack = 0;//counter to look at previous elements for overlap etc.
	int maxcBack = 10;//
	int ndx = 0;
	int tmpStart = 0;
	int tmpStop = 0;
	int maxSP = 0;
	i = j = k = sp = 0;
	BOOLEAN rightShifted = FALSE;
	BOOLEAN leftShifted = FALSE;

	/*******************************************
	 * Start looking for inverted repeats*******
	 *******************************************
	 */
	for (strti = mincrf; strti <= (total_bases - mincrf); strti++) {
		maxSP = min(cspacer,(total_bases-(strti+mincrf)));
		for (sp = 0; sp <= maxSP; sp++) {
			i = strti;
			k = 0;
			j = strti + sp + 1;
			while ((dna[i] == dna3[j]) && (j < (total_bases)) && (i >= 0)
					&& (dna[j] != 'n')) {
				k++;
				j++;
				i--;
			}
			if (k >= mincrf) {
				if ((k <= cut) && (sp > shortSpacer)) {//check for short IR spacers
					continue;
				}
				tmpStart = ((strti - k) + 2);// in ncbi coordinates (+1 to array cords)
				tmpStop = (strti + k + sp + 1);// in ncbi coordinates (+1 to array cords)
				if ((ndx == 0)) {//first one, can't compare current to prev if prev doesn't exist
					rightShifted = FALSE;
					leftShifted = FALSE;
					irep[ndx].start = tmpStart;
					irep[ndx].sub = tmpStop;//min loop boundary set to end by default
					irep[ndx].len = k;
					irep[ndx].loop = sp;
					irep[ndx].num = 1;
					irep[ndx].end = tmpStop;
					irep[ndx].strand = 0;
					ndx++;
				}
				else {//Not first one

					//check for immediate inclusions
					//old within new, new larger looped
					while ((irep[ndx - 1].end <= tmpStop)
							&& (irep[ndx - 1].start >= tmpStart) && irep[ndx
							- 1].len < k && ((ndx - 1) >= 0)) {
						//old within new, new better
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, new better
					while ((irep[ndx - 1].end >= tmpStop)
							&& (irep[ndx - 1].start <= tmpStart) && irep[ndx
							- 1].len < k && ((ndx - 1) >= 0)) {
						ndx--;//replace previous
						rightShifted = FALSE;
						leftShifted = FALSE;

					}
					//old within new, old better
					if ((irep[ndx - 1].end <= tmpStop) && (irep[ndx - 1].start
							>= tmpStart) && irep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					//new within old, old better
					if ((irep[ndx - 1].end >= tmpStop) && (irep[ndx - 1].start
							<= tmpStart) && irep[ndx - 1].len > k) {
						//don't add new
						rightShifted = FALSE;
						leftShifted = FALSE;
					}
					else if ((tmpStop == irep[ndx - 1].end) && (k == irep[ndx
							- 1].len) && (!rightShifted)) {
						leftShifted = TRUE;//to check for alternate shifting
						rightShifted = FALSE;
						ndx--;
						irep[ndx].num = irep[ndx].num + 1;//add one to Permutation count
						irep[ndx].sub = tmpStart;//adjust minimum loop boundary
						ndx++;
					}
					else if ((tmpStart == irep[ndx - 1].start) && (k
							== irep[ndx - 1].len) && (!leftShifted)) {
						rightShifted = TRUE;//to check for alternate shifting
						leftShifted = FALSE;

						//need check as rightShfited grows that it doesn't swallow previous
						if ((irep[ndx - 2].end <= tmpStop)
								&& (irep[ndx - 2].start >= tmpStart)
								&& irep[ndx - 2].len < k) {
							//old within new, new better
							//replace old with current shifting
							irep[ndx - 2] = irep[ndx - 1];
							ndx--;//replace previous
						}

						ndx--;
						irep[ndx].num = (irep[ndx].num + 1);//add one to Permutation count
						irep[ndx].end = tmpStop;//adjust end
						//mrep[ndx].sub = -12;
						ndx++;
					}
					else {//neither shifted, add new repeat
						rightShifted = FALSE;
						leftShifted = FALSE;
						irep[ndx].start = tmpStart;
						irep[ndx].sub = tmpStop;//min loop boundary set to end by default
						irep[ndx].len = k;
						irep[ndx].loop = sp;
						irep[ndx].num = 1;
						irep[ndx].end = tmpStop;
						irep[ndx].strand = 0;
						ndx++;

						for (cBack = 1; cBack <= maxcBack; cBack++) {
							while ((((irep[ndx - (1 + cBack)].end >= irep[ndx
									- 1].end) && (irep[ndx - (1 + cBack)].start
									<= irep[ndx - 1].start)) || ((irep[ndx - (1
									+ cBack)].end <= irep[ndx - 1].end)
									&& (irep[ndx - (1 + cBack)].start
											>= irep[ndx - 1].start)))
									&& ((cBack + 1) <= ndx)) {
								//maximize stem length, then minimize loop length
								if ((irep[ndx - (1 + cBack)].len == irep[ndx
										- 1].len)) {//if stems are equal, keep shortest loop
									//if previous loop is larger, delete it
									if ((irep[ndx - (1 + cBack)].loop
											> irep[ndx - 1].loop)) {
										delIRep((ndx - 1), (ndx - (1 + cBack)));
									}
								}
								//otherwise keep longest stem

								else if ((irep[ndx - (1 + cBack)].len
										< irep[ndx - 1].len)) {//previous stem is shorter = delete it
									delIRep((ndx - 1), (ndx - (1 + cBack)));
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
}/* END of findIR*/
