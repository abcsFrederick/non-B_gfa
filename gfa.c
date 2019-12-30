#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stddef.h>
#include <ctype.h>

#include "gfa.h"

/********************************************************************
 *  Genomic Feature Analyzer; non-B motifs in nucleic acid sequence *
 ********************************************************************/

/*   Structure definition for REP   **************************
 *    typedef struct REP
 *        {
 *        int start;Starting position of repeat
 *		  int loop; Size of Loop (also kv score for Z-DNA)
 *        int pos;  Position of "i", no longer used
 *        int len;  Length of the repeat pattern
 *        int num;  Number of times pattern is repeated (permutations for MR)
 *        int end;  End position of the motif
 *        int ID;   Unique ID for each repeat for the SQL database
 *        char typ; Type of repeat: tandem, mirror, cruciform
 *        int sub;  holds type for str, Islands for G4, MinLoop
 *        	for IR and MR, and remainder for DR
 *        int special; tag for if sequence is Slippped, Cruciform, Triplex etc.
 **************************************************************/

/********************************************************************
 *  declare variables external to avoid the stack overflow problem **
 *  These variables can be up to 300M bytes long!!!	         	   **
 ********************************************************************
 */

//A_Tract atract[16000000];
//short A_Tract_Strt[20000000];
potential_Bent_DNA pAPRs[5 * MAX_REPS + 1];

G_Island gisle[5 * MAX_REPS + 1];
G_Island rcgisle[5 * MAX_REPS + 1]; //reverse comp, actually + strand C islands
const G_Island null_gisle; //defaults to null, used to reset the stuct for each fasta

//potential_G_Quads pGQs[MAX_REPS + 1];
//const potential_G_Quads null_pgq;//defaults to null, used to reset the stuct for each fasta
//potential_G_Quads rcpGQs[MAX_REPS + 1];//reverse comp

//global so we only have to call findGQs and findgislands once
int nGisls;
int nCisls;

char dna[MAX_DNA + 1];
char dna2[MAX_DNA + 1]; //reverse complement strand
char dna3[MAX_DNA + 1]; //complement strand
REP mrep[MAX_REPS + 1]; //mirror
//REP *irep = malloc(2*MAX_REPS * sizeof(REP));
REP irep[MAX_REPS + 1]; //inverted
REP drep[MAX_REPS + 1]; //direct
REP grep[MAX_REPS + 1]; //g-quadraplex
REP zrep[MAX_REPS + 1]; //z-dna
REP srep[MAX_REPS + 1]; //str
REP arep[MAX_REPS + 1]; //a-phased-repeat
const REP null_rep; //defaults to null, used to reset the stuct for each fasta

int main(int argc, char *argv[]) {

	//initialize file name and other strings
	char dna_filename[181];

	char gffout_filenameI[181];
	char gffout_filenameM[181];
	char gffout_filenameD[181];
	char gffout_filenameG[181];
	char gffout_filenameZ[181];
	char gffout_filenameS[181];
	char gffout_filenameA[181];

	char tsvout_filenameI[181];
	char tsvout_filenameM[181];
	char tsvout_filenameD[181];
	char tsvout_filenameG[181];
	char tsvout_filenameZ[181];
	char tsvout_filenameS[181];
	char tsvout_filenameA[181];

	char fasta_title[MAX_FASTA_SIZE + 1];
	char chrom[80];
	char seq_title[80];
	char wget_string[130]; //char array for storing wget command to tell php that prog. is done
	char chmod_stringTSV[130]; //set all to read, set group and owner to write
	char chmod_stringGFF[130]; //set group and owner to write
	strcpy(
			wget_string,
			"wget -S - \"http://nonb.abcc.ncifcrf.gov/modules/nBMSTc/controllers/notify.php?file=");
	strcpy(chmod_stringGFF, "chmod 664 ");
	strcpy(chmod_stringTSV, "chmod 664 ");

	//int w = 0;//fasta title first word counter

	//registered loop counters
	register int i = 0;
	register int j = 0;

	//total bases in input file
	int total_bases = 0;
	int fasta_count = 0;
	int fasta = 0; //main loop counter

	//motif count variables
	int ireps = 0; //n IR, Inverted Repeats	(cruciforms)
	int mreps = 0; //n MR, Mirror Repeats		(triplexes)
	int dreps = 0; //n DR, Direct Repeats		(slipped DNA)
	int greps = 0; //n GQ, G-Quadruplexes		(g-quaduplexes)
	int zreps = 0; //n ZDNA, Z-DNA
	int sreps = 0; //n STR, Short Tandem Repeats
	int areps = 0; //n APR, A-Phased Repeats	(bent DNA)

	//motif precursor counts
	int ngisles = 0; //n g islands
	int rcngisles = 0; //n reverse complement strand g islands
	int nPgreps = 0; //n potential G repeats
	int rcnPgreps = 0; //n reverse complement strand potential G repeats

	//main motif definition default values
	int minGQrep = 3; //no max
	int maxGQspacer = 7; //min =1
	int minMRrep = 10; //no max
	int maxMRspacer = 100; //min =0
	int minIRrep = 6; //no max
	int maxIRspacer = 100; //min =0
	int shortIRcut = 9; //this is the max repeat for for the smaller IRs
	int shortIRspacer = 4; //spacer max for smaller IRs
	int maxDRspacer = 10; //min =0
	int minDRrep = 10;
	int maxDRrep = 300;
	int minATracts = 3; //min number of consecutive a tracks
	int minATractSep = 10; //separation between a tract centers
	int maxATractSep = 11; //separation between a tract centers
	int maxAPRlen = 9;
	int minAPRlen = 3;
	int minZlen = 10; //no max
	int minSTR = 1;
	int maxSTR = 9;
	int minSTRreps = 3;
	int minSTRbp = 10;

	//motif subset definition default values
	int minCruciformRep = 10;
	int maxCruciformSpacer = 4;
	int minTriplexYRpercent = 10;
	int maxTriplexSpacer = 8;
	int maxSlippedSpacer = 0;
	int minKVscore = 33;

	//counters for tsv printing, only want to print header once
	int Ic, Mc, Dc, Zc, Ac, Gc, Sc;
	Ic = Mc = Dc = Zc = Ac = Gc = Sc = 0;

	//booleans to control which motifs/subsets are searched for
	BOOLEAN DO_findSTR, DO_findDR, DO_findMR, DO_findGQ, DO_findIR, DO_findZ,
			DO_findAPR;
	DO_findSTR = DO_findDR = DO_findMR = DO_findGQ = DO_findIR = DO_findZ =
			DO_findAPR = TRUE;

	BOOLEAN DO_Cruciform, DO_Triplex, DO_Slipped, DO_KVzdna;
	DO_Cruciform = DO_Triplex = DO_Slipped = DO_KVzdna = TRUE;
	DO_KVzdna = TRUE;
	BOOLEAN DO_WGET = TRUE; //system call to wget for PHP?
	BOOLEAN DO_CHMOD = FALSE; //system call to change output file permissions, everyone read, owner and group write
	BOOLEAN FATAL = FALSE; //fatal error?
	BOOLEAN ISEQ = FALSE; //input file given?
	//	BOOLEAN GFF = TRUE;//create .gff output file?
	//	BOOLEAN TSV = TRUE;//create .tsv output file?
	BOOLEAN CHROM = FALSE; //chromosome name given as com line argument?
	BOOLEAN KEEP_TIME = TRUE; //for benchmarking etc.

	//time stuff
	time_t startTime;

	FILE *dna_file; //input file name
	//output files
	FILE *gffout_fileI;
	FILE *gffout_fileM;
	FILE *gffout_fileD;
	FILE *gffout_fileG;
	FILE *gffout_fileZ;
	FILE *gffout_fileS;
	FILE *gffout_fileA;

	FILE *tsvout_fileI;
	FILE *tsvout_fileM;
	FILE *tsvout_fileD;
	FILE *tsvout_fileG;
	FILE *tsvout_fileZ;
	FILE *tsvout_fileS;
	FILE *tsvout_fileA;

	/*****************************
	 * Procedure Declarations  ***
	 *****************************
	 */

	void print_usage(char pgmname[]); //prints help text to screen
	//void nulls(char line[], int n); //clears character strings
	void rcdna(int ndna); //computes reverse complement dna
	void cdna(int ndna); //computes complement dna

	//io functions
	int read_fasta(FILE *dna_file, char fasta_title[]);
	int read_mult_fasta(FILE *dna_file, int fasta, char fasta_title[]);
	int get_fasta_count(FILE *dna_file);
	void print_gff_file(FILE *gffout_file, int nreps, char chrom[], char X,
			int total_bases);
	void print_tsv_file(FILE *tsvout_file, int nreps, char chrom[], char X,
			int nFasta, int total_bases);

	//main motif finding functions
	int findIR(int minIRrep, int maxIRspacer, int shortIRcut, int shortIRspacer,
			int total_bases);
	int findMR(int minMRrep, int maxMRspacer, int total_bases);
	int findDR(int minDRrep, int maxDRrep, int maxDRspacer, int total_bases);
	int findZDNA(int minZ, int total_bases);
	int findSTR(int minSTR, int maxSTR, int minSTRbp, int minSTRreps,
			int total_bases);
	int findAPR(int minAPRlen, int maxAPRlen, int minATracts, int total_bases);
	int findGQ(int minGQrep, int maxGQspacer);

	//functions for motif precursor sequences
	int getGislands(int minGQrep, int total_bases);
	int getPuPyruns(int minZ, int total_bases);
	//int process_gislands(int ngisles, int minGQrep, BOOLEAN plus,
	//		int maxGQspacer);
	int getAtracts(int minAT, int maxAT, int total_bases);
//	int process_Atracts(int minAT, int nATs, int total_bases, BOOLEAN plus);

	//motif post processing/filtering functions
	int process_repeatsCentered(int nreps, char X);
	int process_repeatsIncluded(int nreps, char X);
	void is_subset(int nreps, char X, int max_loop, int limit);

	//	char get_sequence(int start, int stop, int rep, char X, int strand){
	/*****************************************
	 ************* START OF PROGRAM  *********
	 *****************************************
	 */

	// if program with no arguments
	if (argc == 1) {
		print_usage(argv[0]);
		exit(1);
	}

	/*******************************
	 * get command line parameters *
	 *******************************
	 */
	for (i = 1; i < argc; i++) {

		//input file name: REQUIRED
		if (strncmp(argv[i], "-seq", 4) == 0) {
			memset((char *) dna_filename, '\0', 180);
			//nulls(dna_filename, 80);
			if (argv[i + 1] != NULL) {
				strncpy(dna_filename, argv[i + 1], strlen(argv[i + 1]));
				fprintf(stderr, "-seq value = %s\n", dna_filename);
				ISEQ = TRUE;
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -seq switch (input file name)\n");
				FATAL = TRUE;
			}
		}
		//Chromosome name, truncated to 10 chars, REQUIRED if .gff output = TRUE
		if (strncmp(argv[i], "-chrom", 6) == 0) {
			memset((char *) chrom, '\0', 10);
			//nulls(chrom, 10);
			if (argv[i + 1] != NULL) {
				j = strlen(argv[i + 1]);
				if (j > 10) {
					j = 10;
				}
				strncpy(chrom, argv[i + 1], j);
				CHROM = TRUE;
				fprintf(stderr, "-chrom value = %s\n", chrom);
			} else {
				fprintf(stderr, "FATAL ERROR: No argument for -chrom switch\n");
				FATAL = TRUE;
			}
		}

		//command line overrides for main motif default values
		if (strncmp(argv[i], "-minGQrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minGQrep);
				fprintf(stderr, "-minGQrep value = %d\n", minGQrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minGQrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxGQspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxGQspacer);
				fprintf(stderr, "-maxGQspacer value = %d\n", maxGQspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxGQspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minMRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minMRrep);
				fprintf(stderr, "-minMRrep value = %d\n", minMRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minMRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minIRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minIRrep);
				fprintf(stderr, "-minIRrep value = %d\n", minIRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minIRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxMRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxMRspacer);
				fprintf(stderr, "-maxMRspacer value = %d\n", maxMRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxMRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxIRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxIRspacer);
				fprintf(stderr, "-maxIRspacer value = %d\n", maxIRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxIRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxDRspacer", 12) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxDRspacer);
				fprintf(stderr, "-maxDRspacer value = %d\n", maxDRspacer);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxDRspacer switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minDRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minDRrep);
				fprintf(stderr, "-minDRrep value = %d\n", minMRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minDRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxDRrep", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxDRrep);
				fprintf(stderr, "-maxDRrep value = %d\n", maxDRrep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxDRrep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minATracts", 11) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minATracts);
				fprintf(stderr, "-minATracts value = %d\n", minATracts);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minATtracts switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minATractSep", 13) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minATractSep);
				fprintf(stderr, "-minATractSep value = %d\n", minATractSep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minATtractSep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxATractSep", 13) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxATractSep);
				fprintf(stderr, "-maxATractSep value = %d\n", maxATractSep);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxATtractSep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxAPRlen", 10) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxAPRlen);
				fprintf(stderr, "-maxAPRlen value = %d\n", maxAPRlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxAPRlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minAPRlen", 10) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minAPRlen);
				fprintf(stderr, "-minAPRlen value = %d\n", minAPRlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minAPRlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minZlen", 8) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minZlen);
				fprintf(stderr, "-minZlen value = %d\n", minZlen);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minZlen switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minSTR", 7) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minSTR);
				fprintf(stderr, "-minSTR value = %d\n", minSTR);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minSTR switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxSTR", 7) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxSTR);
				fprintf(stderr, "-maxSTR value = %d\n", maxSTR);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -maxSTR switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-minSTRbp", 9) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minSTRbp);
				fprintf(stderr, "-minSTRbp value = %d\n", minSTRbp);
			} else {
				fprintf(stderr,
						"FATAL ERROR: No argument for -minSTRbp switch\n");
				FATAL = TRUE;
			}
		}

		//command line overrides for subset motif default values
		if (strncmp(argv[i], "-minCruciformRep", 16) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &minCruciformRep);
				fprintf(stderr, "-minCruciformRep value = %d\n",
						minCruciformRep);
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -maxCruciformRep switch\n");
				FATAL = TRUE;
			}
		}
		if (strncmp(argv[i], "-maxCruciformSpacer", 19) == 0) {
			if (argv[i + 1] != NULL) {
				sscanf(argv[i + 1], "%d", &maxCruciformSpacer);
				fprintf(stderr, "-maxCruciformSpacer value = %d\n",
						maxCruciformSpacer);
			} else {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -maxCruciformSpacer switch\n");
				FATAL = TRUE;
			}
		}

		//command line boolean overrides
		if (strncmp(argv[i], "-skipZ", 6) == 0) {
			DO_findZ = FALSE;
			fprintf(stderr, "-skipZ = Z-DNA motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipDR", 7) == 0) {
			DO_findDR = FALSE;
			fprintf(stderr,
					"-skipDR = Direct Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipIR", 7) == 0) {
			DO_findIR = FALSE;
			fprintf(stderr,
					"-skipIR = Inverted Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipMR", 7) == 0) {
			DO_findMR = FALSE;
			fprintf(stderr,
					"-skipMR = Mirror Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipGQ", 7) == 0) {
			DO_findGQ = FALSE;
			fprintf(stderr,
					"-skipGQ = G-Quadruplex motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipSTR", 8) == 0) {
			DO_findSTR = FALSE;
			fprintf(
					stderr,
					"-skipSTR = Short Tandem Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipAPR", 8) == 0) {
			DO_findAPR = FALSE;
			fprintf(stderr,
					"-skipAPR = A-Phased Repeat motifs will not be found\n");
		}
		if (strncmp(argv[i], "-skipSlipped", 12) == 0) {
			DO_Slipped = FALSE;
			fprintf(
					stderr,
					"-skipSlipped = Slipped subset of Direct Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipTriplex", 12) == 0) {
			DO_Triplex = FALSE;
			fprintf(
					stderr,
					"-skipTriplex = Triples subset of Mirror Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipCruciform", 14) == 0) {
			DO_Cruciform = FALSE;
			fprintf(
					stderr,
					"-skipCruciform = Cruciform subset of Inverted Repeats will not be found\n");
		}
		if (strncmp(argv[i], "-skipKVzdna", 11) == 0) {
			DO_KVzdna = FALSE;
			fprintf(
					stderr,
					"-skipKVzdna = Karen Vasquex subset of Z-DNA will not be found\n");
		}
		if (strncmp(argv[i], "-skipWGET", 9) == 0) {
			DO_WGET = FALSE;
			fprintf(stderr, "-skipWGET = wget PHP call will be skipped\n");
		}
		if (strncmp(argv[i], "-doCHMOD", 10) == 0) {
			DO_CHMOD = TRUE;
			fprintf(stderr,
					"-doCHMOD = output files will get chmod 664 permissions\n");
		}
		//gff
		if (strncmp(argv[i], "-out", 4) == 0) {
			//GFF = TRUE;
			if (argv[i + 1] == NULL) {
				fprintf(
						stderr,
						"FATAL ERROR: No argument for -out switch (output file name prefix)\n");
				FATAL = TRUE;
			} else { //out file name given, create output files
				if (DO_findIR) {
					memset((char *) gffout_filenameI, '\0', 80);
					memset((char *) tsvout_filenameI, '\0', 80);
					strncpy(gffout_filenameI, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameI, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findMR) {
					memset((char *) gffout_filenameM, '\0', 80);
					memset((char *) tsvout_filenameM, '\0', 80);
					strncpy(gffout_filenameM, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameM, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findDR) {
					memset((char *) gffout_filenameD, '\0', 80);
					memset((char *) tsvout_filenameD, '\0', 80);
					strncpy(gffout_filenameD, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameD, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findGQ) {
					memset((char *) gffout_filenameG, '\0', 80);
					memset((char *) tsvout_filenameG, '\0', 80);
					strncpy(gffout_filenameG, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameG, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findZ) {
					memset((char *) gffout_filenameZ, '\0', 80);
					memset((char *) tsvout_filenameZ, '\0', 80);
					strncpy(gffout_filenameZ, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameZ, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findSTR) {
					memset((char *) gffout_filenameS, '\0', 80);
					memset((char *) tsvout_filenameS, '\0', 80);
					strncpy(gffout_filenameS, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameS, argv[i + 1], strlen(argv[i + 1]));
				}
				if (DO_findAPR) {
					memset((char *) gffout_filenameA, '\0', 80);
					memset((char *) tsvout_filenameA, '\0', 80);
					strncpy(gffout_filenameA, argv[i + 1], strlen(argv[i + 1]));
					strncpy(tsvout_filenameA, argv[i + 1], strlen(argv[i + 1]));
				}
				strncat(wget_string, argv[i + 1], strlen(argv[i + 1]));
				strncat(chmod_stringTSV, argv[i + 1], strlen(argv[i + 1]));
				strncat(chmod_stringGFF, argv[i + 1], strlen(argv[i + 1]));
				strncat(wget_string, "\"", 2);
				strncat(chmod_stringTSV, "*.tsv", 5);
				strncat(chmod_stringGFF, "*.gff", 5);

			}

		}

	}

	/*************************************************
	 * check for required command line parameters ****
	 *************************************************
	 */
	if (!ISEQ) {
		fprintf(stderr, " FATAL ERROR: no sequence file given (-seq)\n");
		FATAL = TRUE;
	}
	//if ((GFF && !CHROM) || (TSV && !CHROM)) CHROM_ERROR = TRUE;
	if (!CHROM) {
		fprintf(
				stderr,
				"-chrom not given. Sequence will be named using first word of fasta title line.\n");
	}
	if (FATAL) {
		fprintf(stderr, " FATAL Errors, exiting program\n");
		exit(20);
	}

	//open fasta input file
	if ((dna_file = fopen(dna_filename, "r")) == NULL) {
		fprintf(stderr, "\nERROR - Cannot open dna input file %s \n",
				dna_filename);
		exit(2);
	}

	if (DO_findIR) {
		strncat(gffout_filenameI, "_IR.gff", 7);
		if ((gffout_fileI = fopen(gffout_filenameI, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameI);
			exit(3);
		}
	}
	if (DO_findMR) {
		strncat(gffout_filenameM, "_MR.gff", 7);
		if ((gffout_fileM = fopen(gffout_filenameM, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameM);
			exit(3);
		}
	}
	if (DO_findDR) {
		strncat(gffout_filenameD, "_DR.gff", 7);
		if ((gffout_fileD = fopen(gffout_filenameD, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameD);
			exit(3);
		}
	}
	if (DO_findGQ) {
		strncat(gffout_filenameG, "_GQ.gff", 7);
		if ((gffout_fileG = fopen(gffout_filenameG, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameG);
			exit(3);
		}
	}
	if (DO_findZ) {
		strncat(gffout_filenameZ, "_Z.gff", 6);
		if ((gffout_fileZ = fopen(gffout_filenameZ, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameZ);
			exit(3);
		}
	}
	if (DO_findSTR) {
		strncat(gffout_filenameS, "_STR.gff", 8);
		if ((gffout_fileS = fopen(gffout_filenameS, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameS);
			exit(3);
		}
	}
	if (DO_findAPR) {
		strncat(gffout_filenameA, "_APR.gff", 8);
		if ((gffout_fileA = fopen(gffout_filenameA, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open GFF output file %s \n",
					gffout_filenameA);
			exit(3);
		}
	}

	if (DO_findIR) {
		strncat(tsvout_filenameI, "_IR.tsv", 7);
		if ((tsvout_fileI = fopen(tsvout_filenameI, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameI);
			exit(3);
		}
	}
	if (DO_findMR) {
		strncat(tsvout_filenameM, "_MR.tsv", 7);
		if ((tsvout_fileM = fopen(tsvout_filenameM, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameM);
			exit(3);
		}
	}
	if (DO_findDR) {
		strncat(tsvout_filenameD, "_DR.tsv", 7);
		if ((tsvout_fileD = fopen(tsvout_filenameD, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameD);
			exit(3);
		}
	}
	if (DO_findGQ) {
		strncat(tsvout_filenameG, "_GQ.tsv", 7);
		if ((tsvout_fileG = fopen(tsvout_filenameG, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameG);
			exit(3);
		}
	}
	if (DO_findZ) {
		strncat(tsvout_filenameZ, "_Z.tsv", 6);
		if ((tsvout_fileZ = fopen(tsvout_filenameZ, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameZ);
			exit(3);
		}
	}
	if (DO_findSTR) {
		strncat(tsvout_filenameS, "_STR.tsv", 8);
		if ((tsvout_fileS = fopen(tsvout_filenameS, "w")) == NULL) {
			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameS);
			exit(3);
		}
	}
	if (DO_findAPR) {
		strncat(tsvout_filenameA, "_APR.tsv", 8);
		if ((tsvout_fileA = fopen(tsvout_filenameA, "w")) == NULL) {

			fprintf(stderr, "\nERROR - Cannot open TSV output file %s \n",
					tsvout_filenameA);
			exit(3);
		}
	}

	/*************************************
	 * Initialize and READ dna array   ***
	 * dna array defined global/extern ***
	 * to avoid stack overflow issues  ***
	 *************************************
	 */
	memset((char *) dna, '\0', MAX_DNA);

	fasta_count = get_fasta_count(dna_file);
	fprintf(stderr, "Fasta Sections = %d\n", fasta_count);
	dna_file = fopen(dna_filename, "r");

	/************************************
	 * Master Loop for each FASTA entry *
	 ************************************
	 */
	//for each fasta section
	for (fasta = 1; fasta <= fasta_count; fasta++) {

		//reset dna_file
		fclose(dna_file);
		dna_file = fopen(dna_filename, "r");

		//reset all repeat related counters
		total_bases = 0;

		ireps = 0; //n IR, Inverted Repeats	(cruciforms)
		mreps = 0; //n MR, Mirror Repeats		(triplexes)
		dreps = 0; //n DR, Direct Repeats		(slipped DNA)
		greps = 0; //n GQ, G-Quadruplexes		(g-quaduplexes)
		zreps = 0; //n ZDNA, Z-DNA
		sreps = 0; //n STR, Short Tandem Repeats
		areps = 0; //n APR, A-Phased Repeats	(bent DNA)

		ngisles = 0; //n g islands
		rcngisles = 0; //n reverse complement strand g islands
		nPgreps = 0; //n potential G repeats
		rcnPgreps = 0; //n reverse complement strand potential G repeats

		//read dna
		total_bases = read_mult_fasta(dna_file, fasta, fasta_title);
		if (CHROM) {
			strcpy(seq_title, chrom);
		}

		else {
			int w = 1;
			while (!isspace(fasta_title[w]) && (fasta_title[w] != '\0')) {
				w++;
			}
			memset((char *) seq_title, '\0', 80);
			//nulls(seq_title, 80);
			strncpy(seq_title, &fasta_title[0], w);
		}

		if (CHROM) {
			//append fasta section number
			if (fasta_count > 1) {
				char tmp_str[12];
				memset((char *) tmp_str, '\0', 12);
				//nulls(tmp_str, 12);
				sprintf(tmp_str, "%d", fasta);
				strncat(seq_title, "_", 1);
				strncat(seq_title, tmp_str, 5);
			}
		}

		fprintf(stderr, "Total Bases Read=%d \n", total_bases);
		if (total_bases <= 0) {
			fprintf(stderr, " Empty Sequence!. Exiting!\n");
			fflush(NULL);
			exit(20);
		}
		fflush(NULL);

		/*******************************************
		 * compute complements if necessary ********\
		 *******************************************
		 */

		if (DO_findGQ || DO_findIR) { //get complement dna
			cdna(total_bases);
		}
		if (DO_findGQ || DO_findAPR || DO_findIR) { //get reverse complement dna
			rcdna(total_bases);
		}

		//fprintf(stderr, "Starting Repeat Array Initializations\n");
		/*******************************************
		 * (re)Initialize repeat arrays   **************
		 *******************************************
		 */
		for (i = 0; i < MAX_REPS + 1; i++) {
			irep[i] = null_rep;
			mrep[i] = null_rep;
			drep[i] = null_rep;
			grep[i] = null_rep;
			zrep[i] = null_rep;
			srep[i] = null_rep;
			arep[i] = null_rep;
			gisle[i] = null_gisle;
			//pGQs[i] = null_pgq;
			rcgisle[i] = null_gisle;
			//rcpGQs[i] = null_pgq;
		}
		//fprintf(stderr, "Repeat Arrays Initialized\n");

		/**********************************
		 *** Inverted Repeat Section   ****
		 **********************************
		 */

		if (DO_findIR) {

			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting Inverted Repeat search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting Inverted Repeat search\n");
			}
			ireps = findIR(minIRrep, maxIRspacer, shortIRcut, shortIRspacer,
					total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "IRs found = %d; %s\n", ireps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "IRs found = %d\n", ireps);
			}
			if (DO_Cruciform) {
				fprintf(stderr, "Starting Cruciform search\n");
				is_subset(ireps, 'I', maxCruciformSpacer, minCruciformRep);
				fprintf(stderr, "Done Cruciform search\n\n");
			}
			if (fasta_count == 1) {
				if (ireps > 0) {
					print_gff_file(gffout_fileI, ireps, seq_title, 'I',
							total_bases);
					fclose(gffout_fileI);
					print_tsv_file(tsvout_fileI, ireps, seq_title, 'I', ++Ic,
							total_bases);
					fclose(tsvout_fileI);
				}
			}
		}

		/**********************************
		 ***   G-quad Section   ***********
		 **********************************
		 */ //
		if (DO_findGQ) {

			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting G-Quadruplexe search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting G-Quadruplexe search\n");
			}
			getGislands(minGQrep, total_bases);
			//ngisles = getGislands(minGQrep, total_bases, TRUE); //find g islands in plus strand (gisle[])
			//rcngisles = getGislands(minGQrep, total_bases, FALSE); //find g islands in minus strand (rcgisle[])
			//fprintf(stderr, "ok to here, g/c islands = %d  %d\n",nGisls, nCisls);
			//nPgreps = process_gislands(ngisles, minGQrep, TRUE, maxGQspacer); //plus strand potential g-quads (pGQs[])
			//rcnPgreps = process_gislands(ngisles, minGQrep, FALSE, maxGQspacer); //comp potential g-quads (rcpGQs[])
			//greps = findGQ(nPgreps, rcnPgreps); //combines plus and comp strands, loops through both
			greps = findGQ(minGQrep, maxGQspacer);
			fprintf(stderr, "ok to here, greps = %d \n", greps);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "GQs found = %d; %s\n", greps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "GQs found = %d\n", greps);
			}
			if (fasta_count == 1) {
				if (greps > 0) {
					print_gff_file(gffout_fileG, greps, seq_title, 'G',
							total_bases);
					fclose(gffout_fileG);
					print_tsv_file(tsvout_fileG, greps, seq_title, 'G', ++Gc,
							total_bases);
					fclose(tsvout_fileG);
				}
			}
		}

		/**********************************
		 ***   Mirror Repeat Section   ***********
		 **********************************
		 */ //
		if (DO_findMR) {

			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting Mirror Repeat search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting Mirror Repeat search\n");
			}
			mreps = findMR(minMRrep, maxMRspacer, total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "MRs found = %d; %s\n", mreps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "MRs found = %d\n", mreps);
			}
			if (DO_Triplex) {
				fprintf(stderr, "Starting Triplex search\n");
				is_subset(mreps, 'M', maxTriplexSpacer, minTriplexYRpercent);
				fprintf(stderr, "Done Triplex search\n\n");
			}
			if (fasta_count == 1) {
				if (mreps > 0) {
					print_gff_file(gffout_fileM, mreps, seq_title, 'M',
							total_bases);
					fclose(gffout_fileM);
					print_tsv_file(tsvout_fileM, mreps, seq_title, 'M', ++Mc,
							total_bases);
					fclose(tsvout_fileM);
				}
			}
		}

		/**********************************
		 ***   Direct Repeat Section   ****
		 **********************************
		 */
		if (DO_findDR) {

			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting Direct Repeat search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting Direct Repeat search");
			}
			dreps = findDR(minDRrep, maxDRrep, maxDRspacer, total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "DRs found = %d; %s\n", dreps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "DRs found = %d\n", dreps);
			}
			if (DO_Slipped) {
				fprintf(stderr, "Starting Slipped search\n");
				is_subset(dreps, 'D', maxSlippedSpacer, -999);
				fprintf(stderr, "Done Slipped search\n\n");
			}
			if (fasta_count == 1) {
				if (dreps > 0) {
					print_gff_file(gffout_fileD, dreps, seq_title, 'D',
							total_bases);
					fclose(gffout_fileD);
					print_tsv_file(tsvout_fileD, dreps, seq_title, 'D', ++Dc,
							total_bases);
					fclose(tsvout_fileD);
				}
			}
		}
		/**********************************
		 ***   Z-DNA Section   ************
		 **********************************
		 */
		if (DO_findZ) {
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting Z-DNA search; %s", ctime(&startTime));
			} else {
				fprintf(stderr, "Starting Z-DNA search");
			}
			zreps = findZDNA(minZlen, total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Z-DNA found = %d; %s\n", zreps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "Z-DNA found = %d\n", zreps);
			}
			if (DO_KVzdna) {
				fprintf(stderr, "Starting Karen Vasquex Z-DNA search\n");
				is_subset(zreps, 'Z', -999, minKVscore);
				fprintf(stderr, "Done KV Z-DNA\n\n");
			}
			if (fasta_count == 1) {
				if (zreps > 0) {
					print_gff_file(gffout_fileZ, zreps, seq_title, 'Z',
							total_bases);
					fclose(gffout_fileZ);
					print_tsv_file(tsvout_fileZ, zreps, seq_title, 'Z', ++Zc,
							total_bases);
					fclose(tsvout_fileZ);
				}
			}
		}
		/**********************************
		 ***   Short Tandem Repeat Section   ***********
		 **********************************
		 */
		if (DO_findSTR) {
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting Short Tandem Repeat search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting Short Tandem Repeat search\n");
			}
			sreps = findSTR(minSTR, maxSTR, minSTRbp, minSTRreps, total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "STRs found = %d; %s\n", sreps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "STRs found = %d\n", sreps);
			}
			if (fasta_count == 1) {
				if (sreps > 0) {
					print_gff_file(gffout_fileS, sreps, seq_title, 'S',
							total_bases);
					fclose(gffout_fileS);
					print_tsv_file(tsvout_fileS, sreps, seq_title, 'S', ++Sc,
							total_bases);
					fclose(tsvout_fileS);
				}
			}
		}
		/***********************************
		 ***   A-Phased Repeat Section   ***
		 ***********************************/
		if (DO_findAPR) {
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "Starting A-Phased Repeat search; %s",
						ctime(&startTime));
			} else {
				fprintf(stderr, "Starting A-PHased Repeat search\n");
			}
			areps = findAPR(minAPRlen, maxAPRlen, minATracts, total_bases);
			if (KEEP_TIME) {
				time(&startTime);
				fprintf(stderr, "APRs found = %d; %s\n", areps,
						ctime(&startTime));
			} else {
				fprintf(stderr, "APRs found = %d\n", areps);
			}
			if (fasta_count == 1) {
				if (areps > 0) {
					print_gff_file(gffout_fileA, areps, seq_title, 'A',
							total_bases);
					fclose(gffout_fileA);
					print_tsv_file(tsvout_fileA, areps, seq_title, 'A', ++Ac,
							total_bases);
					fclose(tsvout_fileA);
				}
			}
		}
		/**********************************
		 ***   Output the REPEATS *********
		 **********************************/

		if (fasta_count > 1) { //if multiple fasta, write all at end
			if (DO_findIR) {
				if (ireps > 0)
					print_gff_file(gffout_fileI, ireps, seq_title, 'I',
							total_bases);
			}
			if (DO_findMR) {
				if (mreps > 0)
					print_gff_file(gffout_fileM, mreps, seq_title, 'M',
							total_bases);
			}
			if (DO_findDR) {
				if (dreps > 0)
					print_gff_file(gffout_fileD, dreps, seq_title, 'D',
							total_bases);
			}
			if (DO_findZ) {
				if (zreps > 0)
					print_gff_file(gffout_fileZ, zreps, seq_title, 'Z',
							total_bases);
			}
			if (DO_findGQ) {
				if (greps > 0)
					print_gff_file(gffout_fileG, greps, seq_title, 'G',
							total_bases);
			}
			if (DO_findSTR) {
				if (sreps > 0)
					print_gff_file(gffout_fileS, sreps, seq_title, 'S',
							total_bases);
			}
			if (DO_findAPR) {
				if (areps > 0)
					print_gff_file(gffout_fileA, areps, seq_title, 'A',
							total_bases);
			}
			fprintf(stderr, "GFF write done\n");
			if (DO_findIR) {
				if (ireps > 0)
					print_tsv_file(tsvout_fileI, ireps, seq_title, 'I', ++Ic,
							total_bases);
			}
			if (DO_findMR) {
				if (mreps > 0)
					print_tsv_file(tsvout_fileM, mreps, seq_title, 'M', ++Mc,
							total_bases);
			}
			if (DO_findDR) {
				if (dreps > 0)
					print_tsv_file(tsvout_fileD, dreps, seq_title, 'D', ++Dc,
							total_bases);
			}
			if (DO_findZ) {
				if (zreps > 0)
					print_tsv_file(tsvout_fileZ, zreps, seq_title, 'Z', ++Zc,
							total_bases);
			}
			if (DO_findGQ) {
				if (greps > 0)
					print_tsv_file(tsvout_fileG, greps, seq_title, 'G', ++Gc,
							total_bases);
			}
			if (DO_findSTR) {
				if (sreps > 0)
					print_tsv_file(tsvout_fileS, sreps, seq_title, 'S', ++Sc,
							total_bases);
			}
			if (DO_findAPR) {
				if (areps > 0)
					print_tsv_file(tsvout_fileA, areps, seq_title, 'A', ++Ac,
							total_bases);
			}
			fprintf(stderr, "TSV write done\n");
			//fprintf(stderr, "wget string = %s", wget_string);
		}
	}
	//      	fprintf(stdout, "are we here?\n");
        fclose(dna_file);
	if (fasta_count > 1) {
	  if (DO_findIR) {
		fclose(gffout_fileI);
		fclose(tsvout_fileI);
	  }
	  if (DO_findMR) {
		fclose(gffout_fileM);
		fclose(tsvout_fileM);
	  }
	  if (DO_findDR) {
		fclose(gffout_fileD);
		fclose(tsvout_fileD);
	  }
	  if (DO_findZ) {
		fclose(gffout_fileZ);
		fclose(tsvout_fileZ);
	  }
	  if (DO_findGQ) {
		fclose(gffout_fileG);
		fclose(tsvout_fileG);
	  }
	  if (DO_findSTR) {
		fclose(gffout_fileS);
		fclose(tsvout_fileS);
	  }
	  if (DO_findAPR) {
		fclose(gffout_fileA);
		fclose(tsvout_fileA);

	  }
	}
	if (DO_WGET) {
		system(wget_string);
	}
	if (DO_CHMOD) {
		system(chmod_stringTSV); //set all to read, set group and owner to write
		system(chmod_stringGFF); //set all to read, set group and owner to write
	}

	//free(irep);
	//irep = NULL;
	//	fprintf(stdout, "are we done?\n");
	exit(0);
} /* END of main */
