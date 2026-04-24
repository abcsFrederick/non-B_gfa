#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa.h"

/***********************************************************
 *  procedure to create gff(2) files of repeats           **
 *  ================================================      **
 *	a flat tab-delimited file with 9 columns.  Example:	  **
 *	                                                      **
 * chr6	ABCC	Mirror_Repeat	62094	62118	.	+	. **
 * 	ID=chr6_62094_MR;Repeat=10;Spacer=5;CMP=5R/5Y         **
 *                                                        **
 *      def from http://gmod.org/wiki/GFF2                **
 * reference sequence                                     **
 *   This is the ID of the sequence that is used to       **
 *   establish the coordinate system of the annotation    **
 * The source of the annotation.                          **
 *   This field describes how the annotation was derived. **
 * The annotation method, also known as type.             **
 *   This field describes the type of the annotation,     **
 *   such as "CDS". Together the method and source        **
 *   describe the annotation type.                        **
 * start position                                         **
 *   The start of the annotation relative to the          **
 *   reference sequence.                                  **
 * stop position                                          **
 *   The stop of the annotation relative to the reference **
 *   sequence. Start is always less than or equal to stop.**
 * score                                                  **
 *   For annotations that are associated with a numeric   **
 *   score (for example, a sequence similarity), this     **
 *   field describes the score. The score units are       **
 *   completely unspecified, but for sequence             **
 *   similarities, it is typically percent identity.      **
 *   Annotations that do not have a score can use "."     **
 * strand                                                 **
 *   For those annotations which are strand-specific,     **
 *   this field is the strand on which the annotation     **
 *   resides. It is "+" for the forward strand, "-"       **
 *   for the reverse strand, or "." for annotations       **
 *   that are not stranded.                               **
 * phase                                                  **
 *   For annotations that are linked to proteins, this    **
 *   field describes the phase of the annotation on the   **
 *   codons. It is a number from 0 to 2, or "." for       **
 *   features that have no phase.                         **
 * group                                                  **
 *    GFF provides a simple way of generating annotation  **
 *    hierarchies ("is composed of" relationships) by     **
 *    providing a group field. The group field contains   **
 *    the class and ID of an annotation which is the      **
 *    logical parent of the current one.                  **
 ***********************************************************
 */

//void nulls(char line[], int n);

void print_gff_file(FILE *gff_file, int nreps, char chrom[], char X, int total_bases) {

	BOOLEAN GetSeq = TRUE;//compute ATCG count and print motif sequence

	int max_seq = 100000;
	char Seq[max_seq];
	memset((char *) Seq, '\0', max_seq);
	//nulls(Seq, max_seq);
	int j = 0;

	int start = 0;
	int stop = 0;
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
	 * Section to output the repeats if gff file format       **
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
		fprintf(stderr, " nreps S = %d\n", nreps);
		rep = &srep[0];
		strcpy(rep_type, "Short_Tandem_Repeat");
		strcpy(short_rep_type, "STR");
	}
	fprintf(stderr, " Starting print_GFF_file\n");

	fprintf(gff_file,"##gff-version 3\n");
	for (i = 0; i < nreps; i++) {
		//get composition and sequence
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
			if (rep[i].strand == 0) {
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
				start = total_bases-start-1;
				stop = total_bases-stop-1;
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
		//	k++;
		fprintf(gff_file, "%s", chrom);
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		/* Source */
		fprintf(gff_file, "ABCC");
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		/* Type */
		fprintf(gff_file, "%s", rep_type);
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		/* Start */
		fprintf(gff_file, "%d", rep[i].start);
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		/* stop */
		fprintf(gff_file, "%d", rep[i].end);
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		/* Score */
		fprintf(gff_file, ".");
		/* Tab delimiter */
		fprintf(gff_file, "\t");
		if (X == 'G') {
			/* Strand */
			if (rep[i].strand == 0) {
				fprintf(gff_file, "+");
			}
			else {
				fprintf(gff_file, "-");
			}
			/* Tab delimiter */
			fprintf(gff_file, "\t");
			/* Phase */
			fprintf(gff_file, ".");
			/* Tab delimiter */
			fprintf(gff_file, "\t");
			/* ID */
			fprintf(gff_file, "ID=%s_%d_%d_%s", chrom, rep[i].start,rep[i].end, short_rep_type);
			fprintf(gff_file, ";islands=%d", rep[i].sub);
			fprintf(gff_file, ";runs=%d", rep[i].num);
			fprintf(gff_file, ";max=%d", rep[i].len);

			if (GetSeq) {
				fprintf(gff_file, ";composition=%dA/%dC/%dG/%dT", A, C, G, T);
				fprintf(gff_file, ";sequence=%s", Seq);
			}
		}
		if (X != 'G') {//g-quad
			/* Strand */
			fprintf(gff_file, "+");
			/* Tab delimiter */
			fprintf(gff_file, "\t");
			/* Phase */
			fprintf(gff_file, ".");
			/* Tab delimiter */
			fprintf(gff_file, "\t");
			/* ID */
			fprintf(gff_file, "ID=%s_%d_%d_%s", chrom, rep[i].start,rep[i].end, short_rep_type);
			if (X == 'Z') {//z-dna
				fprintf(gff_file, ";length=%d", rep[i].len);
				fprintf(gff_file, ";score=%d", rep[i].loop);
			}
			else if (X == 'A') {//apr
				fprintf(gff_file, ";tracts=%d", rep[i].num);
			}
			else if (X == 'S') {//str
				fprintf(gff_file, ";length=%d", rep[i].len);
				fprintf(gff_file, ";x%d", rep[i].num);
				fprintf(gff_file, "+%d", rep[i].sub);
				fprintf(gff_file, ";type=%d", rep[i].loop);
			}

			else if (X == 'D') {//dr
				fprintf(gff_file, ";spacer=%d", rep[i].loop);
				fprintf(gff_file, ";repeat=%d", rep[i].len);
				fprintf(gff_file, ";x%d", rep[i].num);
				fprintf(gff_file, "+%d", rep[i].sub);
			}
			else if (X == 'M') {//mr
				fprintf(gff_file, ";spacer=%d", rep[i].loop);
				fprintf(gff_file, ";repeat=%d", rep[i].len);
				fprintf(gff_file, ";perms=%d", rep[i].num);
				fprintf(gff_file, ";minloop=%d", rep[i].sub);
			}
			else {//inverted
				fprintf(gff_file, ";spacer=%d", rep[i].loop);
				fprintf(gff_file, ";repeat=%d", rep[i].len);
				fprintf(gff_file, ";perms=%d", rep[i].num);
				fprintf(gff_file, ";minloop=%d", rep[i].sub);
			}

			if (GetSeq) {
				fprintf(gff_file, ";composition=%dA/%dC/%dG/%dT", A, C, G, T);
				fprintf(gff_file, ";sequence=%s", Seq);
			}
			switch (X) {//subset

				case 'D':
				fprintf(gff_file, ";subset=%d", rep[i].special);//repeats

					break;
				case 'M':
				fprintf(gff_file, ";subset=%d", rep[i].special);//repeats

					break;
				case 'I':
				fprintf(gff_file, ";subset=%d", rep[i].special);//repeats
					break;
				case 'G':
					break;
				case 'Z':
				fprintf(gff_file, ";subset=%d", rep[i].special);//repeats
					break;
				case 'A':
					break;
				case 'S':
					break;
			}
		}
		/* New Line */
		fprintf(gff_file, "\n");
	}
	fprintf(stderr, " Total unique repeats written into GFF file= %d\n", nreps);
	return;
} /* END of print_GFF_file */
