#ifndef GFA_H_
#define GFA_H_

#include <stdlib.h>

/* Signal to gfa.h that FALSE/TRUE are provided by R's Boolean.h */
#ifndef GFA_USING_R
typedef enum { FALSE, TRUE } BOOLEAN;
#else
typedef int BOOLEAN;
#endif

/* Practical caps — used by callers when allocating arrays.
 * Actual sizes use seq_len-based formulas; these are absolute maxima. */
#define MAX_REPS    2500000
#define MAX_DNA   300000000
#define MAXCOL          101
#define MAX_FASTA_SIZE   80
#define MAX_LINE        256

typedef struct REP {
    int start;   /* start position (1-based) */
    int loop;    /* spacer size; KV score for Z-DNA */
    int len;     /* length of repeat unit / G-run size */
    int num;     /* times pattern repeats; permutations for MR */
    int end;     /* end position (1-based) */
    int sub;     /* remainder (DR); min-loop boundary (IR/MR); island count (GQ) */
    int strand;  /* 0 = plus, 1 = minus */
    int special; /* 1 = cruciform/triplex/slipped/KV subset, 0 = no */
} REP;

typedef struct A_Tract {
    int   strt;
    short len;
} A_Tract;

typedef struct potential_Bent_DNA {
    double a_center;
    int    strt;
    int    end;
} potential_Bent_DNA;

typedef struct G_Island {
    int strt;
    int len;
} G_Island;

#ifndef max
#  define max(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef min
#  define min(a,b) (((a)<(b))?(a):(b))
#endif

/* ---- function prototypes ---- */

/* cdna.c / rcdna.c */
void cdna(const char *dna, char *dna3, int ndna);
void rcdna(const char *dna, char *dna2, int ndna);

/* findIR.c */
int findIR(const char *dna, const char *dna3, REP *irep,
           int mincrf, int cspacer, int cut, int shortSpacer,
           int total_bases);

/* findMR.c */
int findMR(const char *dna, REP *mrep,
           int minmir, int mspacer, int total_bases);

/* findDR.c */
int findDR(const char *dna, REP *drep,
           int mindir, int maxdir, int dspacer, int total_bases);

/* findGQ.c */
void getGislands(const char *dna,
                 G_Island *gisle,  int *nGisls,
                 G_Island *rcgisle, int *nCisls,
                 int minGQ, int total_bases);
int  findGQ(G_Island *gisle,  int nGisls,
            G_Island *rcgisle, int nCisls,
            REP *grep, int minGQ, int maxGQspacer);

/* findZDNA.c */
int findZDNA(const char *dna, REP *zrep, int minZ, int total_bases);

/* findSTR.c */
int findSTR(const char *dna, REP *srep,
            int minSTR, int maxSTR, int minSTRlen, int minReps,
            int total_bases);

/* findAPR.c */
int getAtracts(const char *dna, const char *dna2,
               potential_Bent_DNA *pAPRs,
               int minAT, int maxAT, int total_bases);
int findAPR(const char *dna, const char *dna2,
            REP *arep, potential_Bent_DNA *pAPRs,
            int minAPR, int maxAPR, int minATracts, int total_bases);

/* process_repeats.c */
int process_repeatsCentered(REP *rep, int nreps);
int process_repeatsIncluded(REP *rep, int nreps);

/* is_subset.c */
void is_subset(const char *dna, REP *rep, int nreps,
               char X, int max_loop, int limit);

/* nulls.c */
void nulls(char line[], int n);

#endif /* GFA_H_ */
