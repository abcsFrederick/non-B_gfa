#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa.h"

/***********************************************************
 *  procedure to create tab separated motif files         **
 *  ================================================      **
 ***********************************************************
 */

//void nulls(char line[], int n);

void print_tsv_file(FILE *tsv_file, int nreps, char chrom[], char X,
		int nFasta, int total_bases) {

	BOOLEAN GetSeq = TRUE;//compute ATCG count and print motif sequence

	int max_seq = 100000;
	char Seq[max_seq];
	memset((char *) Seq, '\0', max_seq);
	//nulls(Seq, max_seq);
	int j = 0;

	int start;
	int stop;
	int A = 0;
	int C = 0;
	int G = 0;
	int T = 0;

	register int i;
	REP *rep;

	printf("%d\n", nreps);
	printf("%c\n", X);

	char rep_type[100];
	char short_rep_type[6];//very bad things happen if these arn't big enough!!

	/***********************************************************
	 * Section to output the repeats if tsv file format       **
	 ***********************************************************
	 */
	if (X == 'M') {
		rep = &mrep[0];
		strcpy(rep_type, "Mirror_Repeat");
		strcpy(short_rep_type, "MR");

	}
	if (X == 'I') {
		rep = &irep[0];
		strcpy(rep_type, "Inverted_Repeat");
		strcpy(short_rep_type, "IR");
	}
	if (X == 'D') {
		rep = &drep[0];
		strcpy(rep_type, "Direct_Repeat");
		strcpy(short_rep_type, "DR");

	}
	if (X == 'G') {
		rep = &grep[0];
		strcpy(rep_type, "G_Quadruplex_Motif");
		strcpy(short_rep_type, "GQ");
	}
	if (X == 'Z') {
		rep = &zrep[0];
		strcpy(rep_type, "Z_DNA_Motif");
		strcpy(short_rep_type, "ZDNA");
	}
	if (X == 'A') {
		rep = &arep[0];
		strcpy(rep_type, "A_Phased_Repeat");
		strcpy(short_rep_type, "APR");
	}
	if (X == 'S') {
		rep = &srep[0];
		strcpy(rep_type, "Short_Tandem_Repeat");
		strcpy(short_rep_type, "STR");
	}

	fprintf(stderr, " Starting print_TSV_file\n");

	//print text file header
	if (nFasta == 1) {//only if first
		fprintf(
				tsv_file,
				"Sequence_name\tSource\tType\tStart\tStop\tLength\tScore\tStrand\tRepeat\tSpacer\t");

		switch (X) {

			case 'D':
			fprintf(tsv_file, "Repeated");
				break;
			case 'M':
			fprintf(tsv_file, "Permutations");
				break;
			case 'I':
			fprintf(tsv_file, "Permutations");
				break;
			case 'G':
			fprintf(tsv_file, "nIslands/nRuns/maxGQ");
				break;
			case 'Z':
			fprintf(tsv_file, "KVScore");
				break;
			case 'A':
			fprintf(tsv_file, "Tracts");
				break;
			case 'S':
			fprintf(tsv_file, "Repeated");

				break;
		}
		fprintf(tsv_file, "\tSubset");
		if (GetSeq) {
			fprintf(tsv_file, "\tComposition\tSequence");
		}
		fprintf(tsv_file, "\n");
	}

	for (i = 0; i < nreps; i++) {//get composition and sequence
		A = 0;
		C = 0;
		G = 0;
		T = 0;
		if (GetSeq) {
			if (X == 'D' || X == 'M' || X == 'I' || X == 'S') {//give only repeat sequence, not loop
				start = rep[i].start - 1;
				stop = start + rep[i].len;
			}
			else {
				start = rep[i].start - 1;
				stop = rep[i].end - 1;
			}
			if (rep[i].strand == 0) {//plus strand
				strncpy(Seq, &dna[rep[i].start - 1],
						rep[i].end - rep[i].start + 1);
				for (j = start; j < stop; j++) {

					switch (dna[j]) {
						/* just complement in this section */
						case 'a':
						A++;
							break;
						case 'c':
						C++;
							break;
						case 'g':
						G++;
							break;
						case 't':
						T++;
							break;
					}
				}
			}
			else {//negative strand
				start = total_bases - start-1;
				stop = total_bases - stop-1;
				strncpy(Seq, &dna2[stop],
						rep[i].end - rep[i].start + 1);
				for (j = stop; j <= start; j++) {
					switch (dna2[j]) {
						// just complement in this section
						case 'a':
						A++;
							break;
						case 'c':
						C++;
							break;
						case 'g':
						G++;
							break;
						case 't':
						T++;
							break;
					}
				}
			}
			Seq[rep[i].end - rep[i].start + 1] = '\0'; //terminate!
		}

		fprintf(tsv_file, "%s", chrom);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Source */
		fprintf(tsv_file, "ABCC");
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Type */
		fprintf(tsv_file, "%s", rep_type);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Start */
		fprintf(tsv_file, "%d", rep[i].start);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Stop */
		fprintf(tsv_file, "%d", rep[i].end);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Length */
		fprintf(tsv_file, "%d", rep[i].end - rep[i].start + 1);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Score */
		fprintf(tsv_file, "NA");
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Strand */
		if (rep[i].strand == 0) {
			fprintf(tsv_file, "+");
		}
		else {
			fprintf(tsv_file, "-");
		}
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Repeat */
		fprintf(tsv_file, "%d", rep[i].len);
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Spacer */
		if (X == 'G') {
			fprintf(tsv_file, "NA");//G loopS given separately
		}
		else if (X == 'Z') {
			fprintf(tsv_file, "NA");//no z loops
		}
		else {
			fprintf(tsv_file, "%d", rep[i].loop);
		}
		/* Tab delimiter */
		fprintf(tsv_file, "\t");

		/* Details/Score */
		switch (X) {

			case 'D': {
				fprintf(tsv_file, "X%d", rep[i].num);
				fprintf(tsv_file, "+%d", rep[i].sub);
			}
				break;
			case 'M': {
				fprintf(tsv_file, "%d", rep[i].num);//permutations
			}
				break;
			case 'G': {
				fprintf(tsv_file, "%dI/", rep[i].sub);//islands
				fprintf(tsv_file, "%dR/", rep[i].num);//runs
				fprintf(tsv_file, "%dM", rep[i].len);//max
			}
				break;
			case 'I': {
				fprintf(tsv_file, "%d", rep[i].num);//permutations
			}
				break;
			case 'Z': {
				fprintf(tsv_file, "%d", rep[i].loop);//score
			}
				break;
			case 'A': {
				fprintf(tsv_file, "%d", rep[i].num);//tracks
			}
				break;
			case 'S': {
				fprintf(tsv_file, "X%d", rep[i].num);
				fprintf(tsv_file, "+%d", rep[i].sub);
			}
				break;
		}
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		/* Notes */
		switch (X) {
			case 'D':
			fprintf(tsv_file, "%d", rep[i].special);//repeats
				break;
			case 'M':
			fprintf(tsv_file, "%d", rep[i].special);//repeats

				break;
			case 'I':
			fprintf(tsv_file, "%d", rep[i].special);//repeats
				break;
			case 'G':
			fprintf(tsv_file, "NA");//repeats
				break;
			case 'Z':
			fprintf(tsv_file, "%d", rep[i].special);//repeats
				break;
			case 'A':
			fprintf(tsv_file, "NA");//repeats
				break;
			case 'S':
			fprintf(tsv_file, "NA");//repeats
				break;
		}
		/* Tab delimiter */
		fprintf(tsv_file, "\t");
		if (GetSeq) {
			/* Composition */
			fprintf(tsv_file, "%dA/%dC/%dG/%dT", A, C, G, T);
			/* Tab delimiter */
			fprintf(tsv_file, "\t");
			/* Sequence */
			fprintf(tsv_file, "%s", Seq);
		}
		/* New Line */
		fprintf(tsv_file, "\n");
	}
	fprintf(stderr, " Total unique repeats written into TSV file= %d\n", nreps);
	return;
} /* END of print_TSV_file */
