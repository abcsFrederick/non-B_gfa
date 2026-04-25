#include <stdio.h>
void print_usage(char arg0[]) {
	fprintf(stderr, "\n\n");
	fprintf(
			stderr,
			"************************  GFA2    ********************************************\n");
	fprintf(
			stderr,
			"*****************************************************************************\n");
	fprintf(
			stderr,
			"usage:%s -seq <input_fasta_filename> -out <output_file_prefix> [optional_switches]\n",
			arg0);
	fprintf(stderr, " \n");
	fprintf(
			stderr,
			" GFA2 takes in a DNA sequence in fasta format and returns Gene Feature Frormat\n");
	fprintf(stderr,
			" (.gff) and Tab Separated Value (.tsv) files containing \n");
	fprintf(stderr,
			" the location and details of potential non-B DNA forming motifs. \n");
	fprintf(stderr, " \n");
	fprintf(stderr, "Required Switches: \n");
	fprintf(stderr,
			"	-seq <string>; The filename for the input DNA fasta file.\n");
	fprintf(stderr, "	-out <string>; The output filename prefix.\n");
	fprintf(stderr,
			"		Motif abbreviations and file extension are automatically appended.\n");
	fprintf(stderr, " \n");
	fprintf(stderr,
			"Optional Integer Switches:  Each switch is followed by its default value.\n");
	fprintf(stderr, "		All values refer to sequential nucleotides.\n");
	fprintf(stderr,
			"		Note: if an integer switch is given, its associated value is required.\n");
	fprintf(
			stderr,
			"	-minGQrep <3>; The minimum number of consecutive G's to form a G run (no max).\n");
	fprintf(stderr,
			"	-maxGQspacer <7>; The maximum allowed distance between G runs (min of 1).\n");
	fprintf(stderr,
			"	-minMRrep <10>; The minimum length of half of a mirror repeat (no max).\n");
	fprintf(stderr,
			"	-maxMRspacer <100>; The c mirror repeat halves (min = 0).\n");
	fprintf(stderr,
			"	-minIRrep <6>; The minimum length of half of an inverted repeat (no max).\n");
	fprintf(
			stderr,
			"	-maxIRspacer <100>; The maximum allowed distance between inverted repeat halves (min = 0).\n");
	fprintf(
			stderr,
			"	-shortIRcut <9>; The maximum length of half of an inverted repeat for it to be considered \"short\".\n");
	fprintf(
			stderr,
			"	-shortIRspacer <4>; The maximum allowed distance between short inverted repeat halves (min = 0).\n");
	fprintf(stderr,
			"	-minDRrep <10>; The minimum length of half of a direct repeat.\n");
	fprintf(stderr,
			"	-maxDRrep <300>; The maximum length of half of a direct repeat.\n");
	fprintf(
			stderr,
			"	-maxDRspacer <100>; The maximum allowed distance between direct repeat halves (min = 0).\n");
	fprintf(
			stderr,
			"	-minATracts <3>; The minimum number of consecutive A Tracts to form an A-Phased Repeat.\n");
	fprintf(stderr,
			"	-minATractSep <10>; The minimum separation between A Tracts centers.\n");
	fprintf(stderr,
			"	-maxATractSep <11>; The maximum separation between A Tracts centers.\n");
	fprintf(
			stderr,
			"	-maxAPRlen <9>; The maximum number of consecutive As allowed in an A tract.\n");
	fprintf(
			stderr,
			"	-minAPRlen <3>; The minimum number of consecutive As allowed in an A tract.\n");
	fprintf(
			stderr,
			"	-minZlen <10>; The minimum length of Z-DNA alternating purine/pyramadine run (no max).\n");
	fprintf(
			stderr,
			"	-minSTR <1>; The minimum length of repeating element in short tandem repeats.\n");
	fprintf(
			stderr,
			"	-maxSTR <9>; The maximum length of repeating element in short tandem repeats.\n");
	fprintf(
			stderr,
			"	-minSTRbp <8>; The minimum overall length for qualification as a short tandem repeat.\n");
	fprintf(stderr, " \n");
	fprintf(
				stderr,
				"	-minCruciformRep <6>; The minimum repeat length for IR to qualify as cruciform.\n");
	fprintf(
				stderr,
				"	-maxCruciformSpacer <4>; The maximum spacer length for IR to qualify as cruciform.\n");
	fprintf(
				stderr,
				"	-minTriplexYRpercent <10>; The minimum purine/pyramadine percent contend for MR to qualify as triplex.\n");
	fprintf(
				stderr,
				"	-maxTriplexSpacer <8>; The maximum spacer length for MR to qualify as triplex.\n");
	fprintf(
				stderr,
				"	-maxSlippedSpacer <0>; The maximum spacer length for DR to qualify as slipped.\n");


	fprintf(stderr, " \n");
	fprintf(stderr, "Other Optional Switches: (not followed by values)\n");
	fprintf(
			stderr,
			"	-chrom <string>; An identifier for the input sequence, \"chr1\" for example.\n");
	fprintf(stderr,
			"	   If not given, the first word of the fasta title string is used\n");
	fprintf(stderr,
			"	-skipAPR; Do not search for A-Phased Repeats (bent DNA). \n");
	fprintf(stderr, "	-skipSTR; Do not search for Short Tandem Repeats. \n");
	fprintf(stderr,
			"	-skipDR; Do not search for Direct Repeats (slipped DNA). \n");
	fprintf(stderr,
			"	-skipMR; Do not search for Mirror Repeats (triplex DNA). \n");
	fprintf(stderr,
			"	-skipIR; Do not search for Inverted Repeats (cruciform DNA). \n");
	fprintf(stderr, "	-skipGQ; Do not search for G-Quadruplexe motifs. \n");
	fprintf(stderr, "	-skipZ; Do not search for Z DNA motifs. \n");
	fprintf(stderr, "	-skipSlipped; Do not search for slipped subset of DRs. \n");
	fprintf(stderr, "	-skipCruciform; Do not search for cruciform subset of IRs. \n");
	fprintf(stderr, "	-skipTriplex; Do not search for triplex subset of MRs. \n");
	fprintf(
			stderr,
			"	-skipWGET; Do not make wget call to php scripts to signify completion. \n");
	fprintf(stderr,
			"	-doCHMOD; Run a system call to chomd command (664) on output files. \n");
	fprintf(stderr,
			"********************************************************************\n");
	fprintf(stderr, "         EXAMPLE:\n");
	fprintf(
			stderr,
			"%s -seq chr22.fa -out chr22 -chrom chr22 -maxDRrep 20 -minATracts 4 -skipGQ -skipZ \n",
			arg0);
	fprintf(stderr,
			"********************************************************************\n");
	fprintf(stderr, " Author: Duncan E. Donohue, Ph.D.\n");
	fprintf(stderr, " Expansion of work by Jack R. Collins, Ph.D.\n");
	fprintf(stderr,
			" Comments, Bugs, Requests: e-mail donohuede@mail.nih.gov\n");
	fprintf(stderr,
			"********************************************************************\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\n");
	return;
}
