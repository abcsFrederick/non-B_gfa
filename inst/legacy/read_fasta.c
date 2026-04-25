#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "gfa.h"
/**********************************************************
 *
 *      Procedure to read the dna file and put into array
 *
 *
 **********************************************************/
/*
 * Defined in main module file to avoid stack overflow problem
 * can easily be 300Mb
 */
extern char dna[];

int read_fasta(FILE *dna_file, char fasta_title[]) {
	register int i;
	int base, fasta_len;
	char line[MAX_LINE + 1];
	BOOLEAN start;

	/* Declare Procedure */
	//void nulls(char line[], int n);

	fasta_len = MAX_FASTA_SIZE;
	start = TRUE;
	i = 0;

	if ((base = getc(dna_file)) != EOF) {
		if (base != '>') {
			fprintf(stderr, " base read =");
			putc(base, stderr);
			fprintf(stderr, ".\n");
			fprintf(stderr, " FATAL ERROR: Sequence file NOT in FASTA format\n");
			exit(88);
		}
		else
			ungetc(base, dna_file);
	}
	else {
		fprintf(stderr, " End of Sequence file!\n");
		fclose(dna_file);
		return (0);
	}

	while ((base = getc(dna_file)) != EOF) {
		switch (base) {
			case '>': /* if (start) get header information
			 else break and return */
			if (start) {
				if (fgets(line, MAX_LINE, dna_file) != NULL) {
					/* crack some part of header and make title */
					memset((char *) fasta_title, '\0', fasta_len);
					//nulls(fasta_title, fasta_len);
					if (strlen(line) < fasta_len) fasta_len = strlen(line) - 1;
					strncpy(fasta_title, &line[0], fasta_len);
				}
				start = FALSE;
			}
			else {
				ungetc(base, dna_file);
				fprintf(stderr, " INFO: Program read %d bases \n", i);
				return (i);
			}
				break;
			default: /* look for DNA sequence */
			if (isalpha(base)) {
				dna[i] = tolower(base);
				i++;
			}
				break;
		}
	}
	fprintf(stderr, " INFO: Program read %d bases \n", i);
	return (i);
}
