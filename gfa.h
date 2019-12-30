#ifndef GFA_H_
#define GFA_H_
#define MAX_REPS 2500000
//#define MAX_REPSIZE 1500
#define MAX_DNA 300000000
//#define MAX_BINS 50000
#define MAXCOL 101
#define MAX_FASTA_SIZE 80
#define MAX_LINE 256

#undef FALSE
#undef TRUE
typedef enum {
	FALSE, TRUE
} BOOLEAN;

/*   Structure definition for REP   **************************
 *    typedef struct REP
 *        {
 *        int start;Starting position of repeat
 *		  int loop; Size of Loop (also kv score for Z-DNA)
 *        int pos;  Position of "i", internal use only
 *        int len;  Length of the repeat pattern
 *        int num;  Number of times pattern is repeated (permutations for MR)
 *        int end;  End position of the motif
 *        int ID;   Unique ID for each repeat for the SQL database
 *        char typ; Type of repeat: tandem, mirror, cruciform
 *        int sub;  holds type for str, Islands for G4, MinLoop
 *        	for IR and MR, and remainder for DR
 *        int special; tag for if sequence is Slippped, Cruciform, Triplex etc.
 *
 *
 *        }REP;
 **************************************************************/

typedef struct REP {
	int start;
	int loop;
//	int pos;
	int len;
	int num;
	int end;
	int sub;
	int strand; //0 for plus, 1 for minus
	int special;//1 for yes, 0 for no
} REP;



/******************************
 *        A-Tract        *****
 *****************************/
typedef struct A_Tract {
	int strt;
	short len;
} A_Tract;

/******************************
 *  potential Bent DNA      *****
 *****************************/
typedef struct potential_Bent_DNA {
	double a_center;
	int strt;
	int end;
} potential_Bent_DNA;

/******************************
 *        G-Island        *****
 *****************************/
typedef struct G_Island {
	int strt;
	int len;
	//int gap;
} G_Island;

/******************************
 *  potential G-quad      *****
 *****************************/
//typedef struct potential_G_Quads {
//	int island;
//	int islands;
//	int runSize;
//	int pruns1;
//	int pruns2;
//	int pruns3;
//} potential_G_Quads;

//these are global so findGQ can process forward and reverse strand in one pass
extern int nGisls;
extern int nCisls;

extern char dna[MAX_DNA + 1];
extern char dna2[MAX_DNA + 1];//reverse complement DNA
extern char dna3[MAX_DNA + 1];//complement DNA
//extern char mdna[MAX_DNA + 1];

extern G_Island gisle[5*MAX_REPS + 1];
//extern potential_G_Quads pGQs[MAX_REPS + 1];
extern G_Island rcgisle[5*MAX_REPS + 1];
//extern potential_G_Quads rcpGQs[MAX_REPS + 1];

//extern A_Tract	atract[16000000];
//extern short A_Tract_Strt[20000000];

extern potential_Bent_DNA pAPRs[5*MAX_REPS + 1];

//extern REP *irep;



//extern REP *irep = malloc(2*MAX_REPS * sizeof(REP));

//if (irep == NULL) {
//  (void)fprintf(stderr, "ERROR: Malloc failed");
//  (void)exit(EXIT_FAILURE);    /* or return EXIT_FAILURE; */
//}


extern REP irep[MAX_REPS + 1];
extern REP mrep[MAX_REPS + 1];
extern REP drep[MAX_REPS + 1];
extern REP grep[MAX_REPS + 1];
extern REP zrep[MAX_REPS + 1];
extern REP srep[MAX_REPS + 1];
extern REP arep[MAX_REPS + 1];

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

#endif

#endif /* GLOBVARS_H_ */

