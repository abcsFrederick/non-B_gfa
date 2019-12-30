#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

/*******************************************
 * A tract definition:
 *
 *Find all AT stretches of 4 to 9,
 * regardless of A/T order.
 *Find Longest A or AnTn run (NAs).
 *Start at first A.  Restart count at TA steps
 *Find longest run of Ts without preceding A’s  (NTs)
 *NAs-NTs>4?
 *If not, try reverse complement
 *If either are >4, == A tract,
 *return the center bp of the largest ANs
 *
 *******************************************
 */

/******************************
 *  A-Phased Repeat Finder   ****
 *****************************/
int getAtracts(int minAT, int maxAT, int total_bases) {

	int n = 0; //a tract nucleotide counter
	int n_rc = 0;
	int nPATs = 0; //number of processed a tracks

	int maxATlen = 0; //maximum length of A or AnTn pattern
	int maxTlen = 0; //maximum length of Ts not following As
	int Alen = 0; //length of current A stretch
	int Tlen = 0; //length of current T (only) stretch
	int ATlen = 0; //length of current A or AnTn pattern
	int TAlen = 0; //length of current T following A stretch
	int maxATend = 0; //end position of the longest A/AnTn stretch
	int ATend = 0; //end A tract

	//same for reverse complement nucleotides
	int maxATlen_rc = 0; //maximum length of A or AnTn pattern
	int maxTlen_rc = 0; //maximum length of Ts not following As
	int Alen_rc = 0; //length of current A stretch
	int Tlen_rc = 0; //length of current T (only) stretch
	int ATlen_rc = 0; //length of current A or AnTn pattern
	int TAlen_rc = 0; //length of current T following A stretch
	int maxATend_rc = 0; //end position of the longest A/AnTn stretch

	register int t; //a tract counter
	t = 0;

	register int i;
	int nAs = 0;
	//int j = 0;
	i = 0;
	int strt = 0;

	while (i < total_bases) {
		if ((dna[i] == 'a') || (dna[i] == 't')) {
			nAs++;
		}
		else {
			if ((nAs >= minAT) && (nAs <= maxAT)) {

				strt = i - nAs + 1;

				ATend = (strt + nAs); //actually one past end!
				Alen = 0;
				Tlen = 0;
				ATlen = 0;
				maxATlen = 0;
				maxTlen = 0;
				TAlen = 0;
				maxATend = 0;

				Alen_rc = 0;
				Tlen_rc = 0;
				ATlen_rc = 0;
				maxATlen_rc = 0;
				maxTlen_rc = 0;
				TAlen_rc = 0;
				maxATend_rc = 0;

				n_rc = (total_bases - ATend); //convert for reverse comp strand (dna2)
				for (n = strt - 1; n < ATend - 1; n++) { //go through each nucleotide in the A-tract
					n_rc++;
					if (dna[n] == 'a') {
						Tlen = 0;
						TAlen = 0;
						if (dna[n - 1] == 't') {
							Alen = 0;
							ATlen = 0;
						}
						else { // a starts or follows another a
							Alen++;
							ATlen++;
						}
					}
					if (dna2[n_rc] == 'a') {
						Tlen_rc = 0;
						TAlen_rc = 0;
						if (dna2[n_rc - 1] == 't') {
							Alen_rc = 0;
							ATlen_rc = 0;
						}
						else { // a starts or follows another a
							Alen_rc++;
							ATlen_rc++;
						}
					}
					if (dna[n] == 't') {
						if (TAlen < Alen) { //Tn following An stretch
							TAlen++;
							ATlen++;
						}
						else { //T starting or not following An
							Tlen++;
							TAlen = 0;
							ATlen = 0;
							Alen = 0;
						}
					}
					if (dna2[n_rc] == 't') {
						if (TAlen_rc < Alen_rc) { //Tn following An stretch
							TAlen_rc++;
							ATlen_rc++;
						}
						else { //T starting or not following An
							Tlen_rc++;
							TAlen_rc = 0;
							ATlen_rc = 0;
							Alen_rc = 0;
						}
					}
					if (maxATlen < ATlen) {
						maxATlen = ATlen;
						maxATend = n;
					}
					if (maxTlen < Tlen) { //could use max macro...
						maxTlen = Tlen;
					}
					if (maxATlen_rc < ATlen_rc) {
						maxATlen_rc = ATlen_rc;
						maxATend_rc = n_rc;
					}
					if (maxTlen_rc < Tlen_rc) { //could use max macro...
						maxTlen_rc = Tlen_rc;
					}
				}
				//if a tract is valid add to processed A-Phased Repeats (pAPRs)
				if (((maxATlen - maxTlen) >= minAT) || ((maxATlen_rc
						- maxTlen_rc) >= minAT)) {
					pAPRs[nPATs].end = strt + nAs;
					pAPRs[nPATs].strt = strt;
					if ((maxATlen - maxTlen) >= (maxATlen_rc - maxTlen_rc)) {
						pAPRs[nPATs].a_center = ((double) maxATend
								- (((double) maxATlen - 1) / 2)) + 1;
					}
					else {
						pAPRs[nPATs].a_center = total_bases
								- (((double) maxATend_rc
										- (((double) maxATlen_rc - 1) / 2)));
					}
					nPATs++;
				}
			}//end process A-tract loop
			nAs = 0;//reset a/t counter
		}
		i++;
	}
	//	fprintf(stderr, "all a tracts = %d\n", j);
	//	return (j);
	fprintf(stderr, "n potential a tracts = %d\n", nPATs);
	return (nPATs);
}/* END of getAtracts*/



/******************************
 *  A-Phased Repeat Finder   ****
 *****************************/

int findAPR(int minAPR, int maxAPR, int minATracts, int total_bases) {
	register int i;
	//int nATs;
	int nProcessedATs;
	i = 0;
	int nBends;

	nProcessedATs = getAtracts(minAPR, maxAPR, total_bases);


	int tracts = 1;
	double distToNext = 0;
	nBends = 0;
	int ndx = 0;
	for (i = 0; i < nProcessedATs - (minATracts + 1); i++) {
		distToNext = pAPRs[i + 1].a_center - pAPRs[i].a_center;
		if ((distToNext <= 11.1) && (distToNext >= 9.9)) {
			tracts++;
		}
		else {
			if (tracts >= minATracts) {
				arep[ndx].start = pAPRs[(i - tracts) + 1].strt;
				arep[ndx].loop = 0;
				arep[ndx].num = tracts; //number of a-tracts`
				//		arep[ndx].pos = 0;
				arep[ndx].strand = 0;
				arep[ndx].len = tracts;
				arep[ndx].end = pAPRs[i].end - 1;
				ndx++;
			}
			tracts = 1;
		}
	}
	return (ndx);
}
