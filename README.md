# non-B-gfa
gfa programs for Non-B site at NCI/FNLCR

gfa is a Suite of programs developed at NCI-Frederick/Frederick National Lab to find sequences associated with non-B DNA forming motifs

DNA exists in many possible conformations that include the A-DNA, B-DNA, and Z-DNA forms; of these, B-DNA is the most common form found in cells. The DNAs that do not fall into a right-handed Watson-Crick double-helix are known as non-B DNAs and comprise cruciform, triplex, slipped (hairpin) structures, tetraplex (G-quadruplex), left-handed Z-DNA, and others. Several recent publications have provided significant evidence that non-B DNA structures may play a role in DNA instability and mutagenesis, leading to both DNA rearrangements and increased mutational rates, which are hallmark of cancer.

**Website for submitting sequences: https://nonb-abcc.ncifcrf.gov/apps/site/default**

The results from the website are based on the default values of gfa and should match the example output (included in the tar file) when the example (below) is run as shown. 


Please cite: Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools.
Regina Z. Cer, Duncan E. Donohue, Uma S. Mudunuri, Nuri A. Temiz, Michael A. Loss, Nathan J. Starner, Goran N. Halusa, Natalia Volfovsky, Ming Yi, Brian T. Luke, Albino Bacolla, Jack R. Collins and Robert M. Stephens.
Nucl. Acids Res. (2013) 41 (D1): D94-D100. doi: 10.1093/nar/gks955

```
************************  GFA2    ********************************************

usage:./gfa -seq <input_fasta_filename> -out <output_file_prefix> [optional_switches]
*****************************************************************************
 GFA2 takes in a DNA sequence in fasta format and returns Gene Feature Format
 (.gff) and Tab Separated Value (.tsv) files containing the location and details of potential non-B DNA forming motifs. 
 
Required Switches:
	-seq <string>; The filename for the input DNA fasta file.
	-out <string>; The output filename prefix.
	Motif abbreviations and file extension are automatically appended.
 
Optional Integer Switches:  Each switch is followed by its default value.
		All values refer to sequential nucleotides.
		Note: if an integer switch is given, its associated value is required.
	-minGQrep <3>; The minimum number of consecutive G's to form a G run (no max).
	-maxGQspacer <7>; The maximum allowed distance between G runs (min of 1).
	-minMRrep <10>; The minimum length of half of a mirror repeat (no max).
	-maxMRspacer <100>; The c mirror repeat halves (min = 0).
	-minIRrep <6>; The minimum length of half of an inverted repeat (no max).
	-maxIRspacer <100>; The maximum allowed distance between inverted repeat halves (min = 0).
	-shortIRcut <9>; The maximum length of half of an inverted repeat for it to be considered "short".
	-shortIRspacer <4>; The maximum allowed distance between short inverted repeat halves (min = 0).
	-minDRrep <10>; The minimum length of half of a direct repeat.
	-maxDRrep <300>; The maximum length of half of a direct repeat.
	-maxDRspacer <100>; The maximum allowed distance between direct repeat halves (min = 0).
	-minATracts <3>; The minimum number of consecutive A Tracts to form an A-Phased Repeat.
	-minATractSep <10>; The minimum separation between A Tracts centers.
	-maxATractSep <11>; The maximum separation between A Tracts centers.
	-maxAPRlen <9>; The maximum number of consecutive As allowed in an A tract.
	-minAPRlen <3>; The minimum number of consecutive As allowed in an A tract.
	-minZlen <10>; The minimum length of Z-DNA alternating purine/pyramadine run (no max).
	-minSTR <1>; The minimum length of repeating element in short tandem repeats.
	-maxSTR <9>; The maximum length of repeating element in short tandem repeats.
	-minSTRbp <8>; The minimum overall length for qualification as a short tandem repeat.
	-minCruciformRep <6>; The minimum repeat length for IR to qualify as cruciform.
	-maxCruciformSpacer <4>; The maximum spacer length for IR to qualify as cruciform.
	-minTriplexYRpercent <10>; The minimum purine/pyramadine percent contend for MR to qualify as triplex.
	-maxTriplexSpacer <8>; The maximum spacer length for MR to qualify as triplex.
	-maxSlippedSpacer <0>; The maximum spacer length for DR to qualify as slipped.
 
Other Optional Switches: (not followed by values)
	-chrom <string>; An identifier for the input sequence, "chr1" for example.
	   If not given, the first word of the fasta title string is used
	-skipAPR; Do not search for A-Phased Repeats (bent DNA). 
	-skipSTR; Do not search for Short Tandem Repeats. 
	-skipDR; Do not search for Direct Repeats (slipped DNA). 
	-skipMR; Do not search for Mirror Repeats (triplex DNA). 
	-skipIR; Do not search for Inverted Repeats (cruciform DNA). 
	-skipGQ; Do not search for G-Quadruplexe motifs. 
	-skipZ; Do not search for Z DNA motifs. 
	-skipSlipped; Do not search for slipped subset of DRs. 
	-skipCruciform; Do not search for cruciform subset of IRs. 
	-skipTriplex; Do not search for triplex subset of MRs. 
	-skipWGET; Do not make wget call to php scripts to signify completion. 
	-doCHMOD; Run a system call to chomd command (664) on output files. 
********************************************************************
         EXAMPLE:
./gfa -skipWGET -seq gfa_test.fasta -out gfa_test
	The input sequence file is gfa_test.fasta
	There should be 14 output files (included in test_files.tar for comparison) using the default values
********************************************************************
 Author: Duncan E. Donohue, Ph.D.
 Expansion of work by Jack R. Collins, Ph.D.
```
