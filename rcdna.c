#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

//extern char dna[];
//extern char dna2[];

void rcdna(int ndna) {

	register int i, k;

	/*******************************
	 Form the reverse complements **
	 ******************************/
	for (i = 0; i < ndna; i++) {
		k = ndna - i - 1;
		switch (dna[i]) {
			/* just complement in this section */
			case 'a':
			dna2[k] = 't';
				break;
			case 'c':
			dna2[k] = 'g';
				break;
			case 'g':
			dna2[k] = 'c';
				break;
			case 't':
			dna2[k] = 'a';
				break;
			case 'n':
			dna2[k] = 'n';
				break;
		}
	}
	/****************************
	 End of reverse complements **
	 ****************************/
	fprintf(stderr, "Reverse Complement Finished:%d\n", i);
	return;
} /* END */
