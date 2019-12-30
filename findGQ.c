#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "gfa.h"

/*******************************************
 * Start looking for G_Islands *************
 *******************************************
 */
void getGislands(int minGQ, int total_bases) {
	register int i;
	int ngs;//number of consecutive g's
	int ncs;//number of consecutive c's

	nGisls = 0;//global g island index/count
	nCisls = 0;//global c island index/count
	i = 0;
	ngs = 0;
	ncs = 0;
	while (i <= total_bases) {
		//consecutive g's
		if (dna[i] == 'g') {
			ngs++;
		}
		else {
			if (ngs >= minGQ) {
				gisle[nGisls].strt = i - ngs + 1; //+1 because dna array starts at 0?
				gisle[nGisls].len = ngs;
				nGisls++;
			}
			ngs = 0;
		}
		//consecutive c's
		if (dna[i] == 'c') {
			ncs++;
		}
		else {
			if (ncs >= minGQ) {
				rcgisle[nCisls].strt = i - ncs + 1; //+1 because dna array starts at 0?
				rcgisle[nCisls].len = ncs;
				nCisls++;
			}
			ncs = 0;
		}
		i++;
	}
}/* END of getGislands*/

/*******************************************
 * Process G_Islands to look for ***********
 * G-quadruplex forming motifs *************
 *******************************************
 */
int findGQ(int minGQ, int maxGQspacer) {
	int ndx = 0;
	int nIls;
	nIls = 0;
	G_Island *islands;
	register int i, i2;
	int npos, j, k, m;
	i = j = k = m = i2 = 0;
	int nposMax, maxGQ;
	nposMax = maxGQ = 0;
	npos = 0;
	int conIls;//consecutive islands
	int strand;
	for (strand = 0; strand < 2; strand++) { //do once for each strand
		if (strand == 0) {//plus strand
			nIls = nGisls;
			islands = &gisle[0];
		}
		if (strand == 1) {//rc strand
			nIls = nCisls;
			islands = &rcgisle[0];
		}
		for (i = 0; i < nIls; i++) {
			conIls = 1;
			npos = (int) (floor((islands[i].len + 1) / (minGQ + 1))); //how many runs of min size will fit in island i?
			i2 = i + 1;
			while (((islands[i2].strt - (islands[i2 - 1].strt
					+ islands[i2 - 1].len)) <= maxGQspacer) && (i2 < nIls)) {//next island is close enough
				conIls++;
				npos += (int) (floor((islands[i2].len + 1) / (minGQ + 1)));
				i2++;
			}
			if (npos >= 4) {
				//seperate loop to find largest possible run
				maxGQ = minGQ;
				for (j = i; j<i2; j++) {
					for (k = islands[j].len; k>maxGQ; k--) {//count down from largest possible in island
						nposMax = (int) (floor((islands[j].len + 1) / (k + 1)));
						for (m = j+1; m<i2; m++) {//through rest of islands in current GQ motif
							nposMax+= (int) (floor((islands[m].len + 1) / (k + 1)));
							if (nposMax>=4) {
								maxGQ = k;
								break;
							}
							if ((int) (floor((islands[m].len + 1) / (k + 1)))==0) {
								if(islands[m+1].strt>(islands[m-1].strt + islands[m-1].len + maxGQspacer)) {//next island is not close enough
									break;
								}
							}
						}
					}
				}
				grep[ndx].start = islands[i].strt;
				grep[ndx].num = npos; //number of min size GQ runs that fit in all islands
				grep[ndx].sub = conIls;//number of islands
				grep[ndx].len = maxGQ;//max GQ possible, separate loop below to raise this if possible
				grep[ndx].end = (islands[i2 - 1].strt + islands[i2 - 1].len)
						- 1;
				grep[ndx].strand = strand;
				ndx++;
			}
			i = i + conIls - 1;//can skip rest of islands in this GQ
		}
	}
	return (ndx);
}//end of findGQ
