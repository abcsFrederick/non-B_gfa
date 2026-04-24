/* Must define GFA_USING_R before including gfa.h so it skips the
 * BOOLEAN enum (R.h already defines FALSE/TRUE as enum constants in C23). */
#include <R.h>
#include <Rinternals.h>
#define GFA_USING_R
#include "gfa.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

/* Lowercase an R character string into a newly malloc'd C string.
 * Caller must free() the result. */
static char *r_seq_to_lower(SEXP r_seq, int *seq_len_out) {
    const char *src = CHAR(STRING_ELT(r_seq, 0));
    int n = (int)strlen(src);
    char *dst = (char *)malloc(n + 1);
    if (!dst) error("nonbgfa: out of memory allocating DNA buffer");
    int i;
    for (i = 0; i < n; i++) dst[i] = (char)tolower((unsigned char)src[i]);
    dst[n] = '\0';
    *seq_len_out = n;
    return dst;
}

/* Convert a REP array to a named R list (one vector per field).
 * The R layer wraps this in as.data.frame(). */
static SEXP rep_to_list(REP *rep, int n) {
    int i;
    SEXP result    = PROTECT(allocVector(VECSXP,  8));
    SEXP names     = PROTECT(allocVector(STRSXP,  8));
    SEXP s_start   = PROTECT(allocVector(INTSXP,  n));
    SEXP s_end     = PROTECT(allocVector(INTSXP,  n));
    SEXP s_strand  = PROTECT(allocVector(STRSXP,  n));
    SEXP s_length  = PROTECT(allocVector(INTSXP,  n));
    SEXP s_spacer  = PROTECT(allocVector(INTSXP,  n));
    SEXP s_num     = PROTECT(allocVector(INTSXP,  n));
    SEXP s_rem     = PROTECT(allocVector(INTSXP,  n));
    SEXP s_subset  = PROTECT(allocVector(LGLSXP,  n));

    for (i = 0; i < n; i++) {
        INTEGER(s_start)[i]  = rep[i].start;
        INTEGER(s_end)[i]    = rep[i].end;
        SET_STRING_ELT(s_strand, i, mkChar(rep[i].strand == 0 ? "+" : "-"));
        INTEGER(s_length)[i] = rep[i].len;
        INTEGER(s_spacer)[i] = rep[i].loop;
        INTEGER(s_num)[i]    = rep[i].num;
        INTEGER(s_rem)[i]    = rep[i].sub;
        LOGICAL(s_subset)[i] = rep[i].special;
    }

    SET_VECTOR_ELT(result, 0, s_start);
    SET_VECTOR_ELT(result, 1, s_end);
    SET_VECTOR_ELT(result, 2, s_strand);
    SET_VECTOR_ELT(result, 3, s_length);
    SET_VECTOR_ELT(result, 4, s_spacer);
    SET_VECTOR_ELT(result, 5, s_num);
    SET_VECTOR_ELT(result, 6, s_rem);
    SET_VECTOR_ELT(result, 7, s_subset);

    SET_STRING_ELT(names, 0, mkChar("start"));
    SET_STRING_ELT(names, 1, mkChar("end"));
    SET_STRING_ELT(names, 2, mkChar("strand"));
    SET_STRING_ELT(names, 3, mkChar("length"));
    SET_STRING_ELT(names, 4, mkChar("spacer"));
    SET_STRING_ELT(names, 5, mkChar("num_repeats"));
    SET_STRING_ELT(names, 6, mkChar("remainder"));
    SET_STRING_ELT(names, 7, mkChar("subset"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(10);
    return result;
}

/* -----------------------------------------------------------------------
 * .Call() entry points — one per motif type
 * ----------------------------------------------------------------------- */

SEXP gfa_find_ir(SEXP r_seq,
                 SEXP r_minIRrep, SEXP r_maxIRspacer,
                 SEXP r_shortIRcut, SEXP r_shortIRspacer,
                 SEXP r_minCruciformRep, SEXP r_maxCruciformSpacer) {
    int seq_len = 0;
    char *dna = r_seq_to_lower(r_seq, &seq_len);

    int cap = max(1024, seq_len / 12);
    char *dna3   = (char *)malloc(seq_len + 1);
    REP  *irep   = (REP  *)calloc(cap, sizeof(REP));
    if (!dna3 || !irep) { free(dna); free(dna3); free(irep); error("nonbgfa: out of memory"); }

    cdna(dna, dna3, seq_len);

    int n = findIR(dna, dna3, irep,
                   asInteger(r_minIRrep), asInteger(r_maxIRspacer),
                   asInteger(r_shortIRcut), asInteger(r_shortIRspacer),
                   seq_len);

    is_subset(dna, irep, n, 'I',
              asInteger(r_maxCruciformSpacer), asInteger(r_minCruciformRep));

    SEXP result = PROTECT(rep_to_list(irep, n));
    free(dna); free(dna3); free(irep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_mr(SEXP r_seq,
                 SEXP r_minMRrep, SEXP r_maxMRspacer,
                 SEXP r_minTriplexYRpercent, SEXP r_maxTriplexSpacer) {
    int seq_len = 0;
    char *dna = r_seq_to_lower(r_seq, &seq_len);

    int cap = max(1024, seq_len / 20);
    REP *mrep = (REP *)calloc(cap, sizeof(REP));
    if (!mrep) { free(dna); error("nonbgfa: out of memory"); }

    int n = findMR(dna, mrep,
                   asInteger(r_minMRrep), asInteger(r_maxMRspacer),
                   seq_len);

    is_subset(dna, mrep, n, 'M',
              asInteger(r_maxTriplexSpacer), asInteger(r_minTriplexYRpercent));

    SEXP result = PROTECT(rep_to_list(mrep, n));
    free(dna); free(mrep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_dr(SEXP r_seq,
                 SEXP r_minDRrep, SEXP r_maxDRrep, SEXP r_maxDRspacer,
                 SEXP r_maxSlippedSpacer) {
    int seq_len = 0;
    char *dna = r_seq_to_lower(r_seq, &seq_len);

    int cap = max(1024, seq_len / 20);
    REP *drep = (REP *)calloc(cap, sizeof(REP));
    if (!drep) { free(dna); error("nonbgfa: out of memory"); }

    int n = findDR(dna, drep,
                   asInteger(r_minDRrep), asInteger(r_maxDRrep),
                   asInteger(r_maxDRspacer), seq_len);

    is_subset(dna, drep, n, 'D', asInteger(r_maxSlippedSpacer), -999);

    SEXP result = PROTECT(rep_to_list(drep, n));
    free(dna); free(drep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_gq(SEXP r_seq, SEXP r_minGQrep, SEXP r_maxGQspacer) {
    int seq_len = 0;
    int minGQ   = asInteger(r_minGQrep);
    char *dna   = r_seq_to_lower(r_seq, &seq_len);

    int island_cap = max(1024, seq_len / minGQ + 1);
    int rep_cap    = max(1024, seq_len / 12);
    int nGisls = 0, nCisls = 0;

    G_Island *gisle   = (G_Island *)calloc(island_cap, sizeof(G_Island));
    G_Island *rcgisle = (G_Island *)calloc(island_cap, sizeof(G_Island));
    REP      *grep    = (REP      *)calloc(rep_cap,    sizeof(REP));
    if (!gisle || !rcgisle || !grep) {
        free(dna); free(gisle); free(rcgisle); free(grep);
        error("nonbgfa: out of memory");
    }

    getGislands(dna, gisle, &nGisls, rcgisle, &nCisls, minGQ, seq_len);
    int n = findGQ(gisle, nGisls, rcgisle, nCisls, grep,
                   minGQ, asInteger(r_maxGQspacer));

    SEXP result = PROTECT(rep_to_list(grep, n));
    free(dna); free(gisle); free(rcgisle); free(grep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_zdna(SEXP r_seq, SEXP r_minZlen, SEXP r_minKVscore) {
    int seq_len = 0;
    char *dna = r_seq_to_lower(r_seq, &seq_len);

    int cap = max(1024, seq_len / 10);
    REP *zrep = (REP *)calloc(cap, sizeof(REP));
    if (!zrep) { free(dna); error("nonbgfa: out of memory"); }

    int n = findZDNA(dna, zrep, asInteger(r_minZlen), seq_len);
    is_subset(dna, zrep, n, 'Z', -999, asInteger(r_minKVscore));

    SEXP result = PROTECT(rep_to_list(zrep, n));
    free(dna); free(zrep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_str(SEXP r_seq,
                  SEXP r_minSTR, SEXP r_maxSTR,
                  SEXP r_minSTRbp, SEXP r_minSTRreps) {
    int seq_len = 0;
    char *dna = r_seq_to_lower(r_seq, &seq_len);

    int cap = max(1024, seq_len / 8);
    REP *srep = (REP *)calloc(cap, sizeof(REP));
    if (!srep) { free(dna); error("nonbgfa: out of memory"); }

    int n = findSTR(dna, srep,
                    asInteger(r_minSTR), asInteger(r_maxSTR),
                    asInteger(r_minSTRbp), asInteger(r_minSTRreps),
                    seq_len);

    SEXP result = PROTECT(rep_to_list(srep, n));
    free(dna); free(srep);
    UNPROTECT(1);
    return result;
}

SEXP gfa_find_apr(SEXP r_seq,
                  SEXP r_minAPRlen, SEXP r_maxAPRlen,
                  SEXP r_minATracts) {
    int seq_len = 0;
    int minAPR  = asInteger(r_minAPRlen);
    char *dna   = r_seq_to_lower(r_seq, &seq_len);

    int apr_cap  = max(1024, seq_len / 30);
    int pAPR_cap = max(1024, seq_len / minAPR + 1);
    char              *dna2  = (char *)malloc(seq_len + 1);
    REP               *arep  = (REP  *)calloc(apr_cap,  sizeof(REP));
    potential_Bent_DNA *pAPRs = (potential_Bent_DNA *)calloc(pAPR_cap,
                                                    sizeof(potential_Bent_DNA));
    if (!dna2 || !arep || !pAPRs) {
        free(dna); free(dna2); free(arep); free(pAPRs);
        error("nonbgfa: out of memory");
    }

    rcdna(dna, dna2, seq_len);

    int n = findAPR(dna, dna2, arep, pAPRs,
                    minAPR, asInteger(r_maxAPRlen),
                    asInteger(r_minATracts), seq_len);

    SEXP result = PROTECT(rep_to_list(arep, n));
    free(dna); free(dna2); free(arep); free(pAPRs);
    UNPROTECT(1);
    return result;
}
