

# nonbgfa

<!-- badges: start -->

[![R-CMD-check](https://github.com/BIFX547-26/non-b-gfa-johnsonra/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/BIFX547-26/non-b-gfa-johnsonra/actions/workflows/R-CMD-check.yml)
<!-- badges: end -->

**nonbgfa** is an R package for finding non-B DNA-forming motifs in
genomic sequences. It wraps the [GFA
suite](https://nonb-abcc.ncifcrf.gov/apps/site/default) developed at
NCI-Frederick / Frederick National Laboratory for Cancer Research,
exposing the C algorithms via a `.Call()` interface so results are
returned directly as R `data.frame` or `GenomicRanges` objects.

Non-B DNA structures deviate from the canonical right-handed
Watson–Crick double helix and are implicated in genomic instability, DNA
recombination hotspots, and elevated mutation rates in cancer.

## Motif types

| Function      | Motif                | Structural form    |
|---------------|----------------------|--------------------|
| `find_ir()`   | Inverted repeats     | Cruciform DNA      |
| `find_mr()`   | Mirror repeats       | Triplex (H-DNA)    |
| `find_dr()`   | Direct repeats       | Slipped-strand DNA |
| `find_gq()`   | G-quadruplexes       | G4 / tetraplex DNA |
| `find_zdna()` | Z-DNA                | Left-handed helix  |
| `find_str()`  | Short tandem repeats | Microsatellites    |
| `find_apr()`  | A-phased repeats     | Bent DNA           |

## Installation

``` r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("abcsFrederick/non-b-gfa", ref = 'Rpkg')
```

The package requires a C compiler (provided by
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows,
Xcode command-line tools on macOS, or `build-essential` on Linux).

## Quick start

``` r
library(nonbgfa)

# Use the bundled example sequence
fasta <- system.file("extdata", "gfa_test.fasta", package = "nonbgfa")

# Run all seven finders at once
results <- find_nonb(fasta)
sapply(results, nrow)
```

     IR  MR  DR  GQ   Z STR APR 
     14   5   4   7   4   7   1 

``` r
# Or call individual finders
ir  <- find_ir(fasta)
gq  <- find_gq(fasta)
head(ir)
```

|        | seq_name | start |  end | strand | length | spacer | num_repeats | remainder | subset |
|:-------|:---------|------:|-----:|:-------|-------:|-------:|------------:|----------:|:-------|
| Test.1 | seq1     |    16 |   29 | \+     |      6 |      2 |           1 |        29 | TRUE   |
| Test.2 | seq1     |   108 |  121 | \+     |      6 |      2 |           1 |       121 | TRUE   |
| Test.3 | seq1     |   542 |  555 | \+     |      6 |      2 |           1 |       555 | TRUE   |
| Test.4 | seq1     |   813 |  826 | \+     |      6 |      2 |           1 |       826 | TRUE   |
| Test.5 | seq1     |  2251 | 2267 | \+     |      7 |      3 |           1 |      2267 | TRUE   |
| Test.6 | seq1     |  2370 | 2391 | \+     |     10 |      2 |           1 |      2391 | TRUE   |

## Output columns

Every finder returns a `data.frame` with these columns:

| Column | Description |
|----|----|
| `seq_name` | Sequence identifier from the FASTA `>` header |
| `start` | Start position (1-based, inclusive) |
| `end` | End position (1-based, inclusive) |
| `strand` | `"+"` (sense) or `"-"` (antisense) |
| `length` | Repeat-unit length (bp); G-run size for GQ |
| `spacer` | Spacer between repeat halves; KV score for Z-DNA |
| `num_repeats` | Number of times the unit repeats; permutations for MR |
| `remainder` | Partial repeat (DR); min-loop boundary (IR/MR); island count (GQ) |
| `subset` | `TRUE` if locus qualifies as cruciform/triplex/slipped/high-KV |

Pass `format = "GRanges"` to any finder for a
[`GenomicRanges::GRanges`](https://bioconductor.org/packages/GenomicRanges)
object (requires the **GenomicRanges** Bioconductor package).

``` r
find_ir(fasta, format = 'GRanges') |>
    head()
```

    GRanges object with 6 ranges and 5 metadata columns:
          seqnames    ranges strand |    length    spacer num_repeats remainder
             <Rle> <IRanges>  <Rle> | <integer> <integer>   <integer> <integer>
      [1]     seq1     16-29      + |         6         2           1        29
      [2]     seq1   108-121      + |         6         2           1       121
      [3]     seq1   542-555      + |         6         2           1       555
      [4]     seq1   813-826      + |         6         2           1       826
      [5]     seq1 2251-2267      + |         7         3           1      2267
      [6]     seq1 2370-2391      + |        10         2           1      2391
             subset
          <logical>
      [1]      TRUE
      [2]      TRUE
      [3]      TRUE
      [4]      TRUE
      [5]      TRUE
      [6]      TRUE
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Parameter defaults

All parameters match the original `gfa` CLI defaults and the Non-B DB
website.

| Motif | Parameter | Default | Meaning |
|----|----|----|----|
| IR | `minIRrep` | 6 | Min arm length (bp) |
| IR | `maxIRspacer` | 100 | Max spacer (bp) |
| IR | `shortIRcut` | 9 | Arms ≤ this length require tight spacer |
| IR | `shortIRspacer` | 4 | Spacer limit for short IRs |
| MR | `minMRrep` | 10 | Min half-length (bp) |
| MR | `maxMRspacer` | 100 | Max spacer (bp) |
| DR | `minDRrep` / `maxDRrep` | 10 / 300 | Repeat unit size range (bp) |
| DR | `maxDRspacer` | 10 | Max spacer (bp) |
| GQ | `minGQrep` | 3 | Min consecutive G’s per run |
| GQ | `maxGQspacer` | 7 | Max spacer between G-runs (bp) |
| Z-DNA | `minZlen` | 10 | Min alternating pur/pyr run (bp) |
| STR | `minSTR` / `maxSTR` | 1 / 9 | Repeat unit size range (bp) |
| STR | `minSTRbp` | 10 | Min total locus length (bp) |
| APR | `minAPRlen` / `maxAPRlen` | 3 / 9 | A-tract length range (bp) |
| APR | `minATracts` | 3 | Min A-tracts per repeat |

See `?find_nonb` and individual function help pages for the full
parameter list.

## Multi-sequence input

Functions will accept either a single named character vector or a
multi-sequence FASTA file:

``` r
# input from fasta file
fasta_ir <- read_fasta(fasta) |>
    find_ir()

# Inline sequences
inline <- c(my_gene = "atcgatcgatcgatcgatcg") |>
    find_ir()
```

## Citation

Please cite the original GFA tool:

> Cer RZ, Donohue DE, Mudunuri US, Temiz NA, Loss MA, Starner NJ, Halusa
> GN, Volfovsky N, Yi M, Luke BT, Bacolla A, Collins JR, Stephens RM.
> (2013) Non-B DB v2.0: a database of predicted non-B DNA-forming motifs
> and its associated tools. *Nucleic Acids Research*, 41(D1):D94–D100.
> <https://doi.org/10.1093/nar/gks955>

## Original C tool

The original `gfa` command-line tool source and Makefile are preserved
in `inst/legacy/` for reference. To build the standalone binary:

``` sh
cd inst/legacy
make
./gfa -skipWGET -seq ../../inst/extdata/gfa_test.fasta -out gfa_test
```
