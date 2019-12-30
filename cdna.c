#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "gfa.h"

//extern char dna[];
//extern char dna2[];

void cdna(int ndna) {

	register int i;

	/*******************************
	 Form the reverse complements **
	 ******************************/
	for (i = 0; i < ndna; i++) {
		switch (dna[i]) {
			/* just complement in this section */
			case 'a':
			dna3[i] = 't';
				break;
			case 'c':
			dna3[i] = 'g';
				break;
			case 'g':
			dna3[i] = 'c';
				break;
			case 't':
			dna3[i] = 'a';
				break;
			case 'n':
			dna3[i] = 'n';
				break;
		}
	}
	/****************************
	 End of reverse complements **
	 ****************************/
	fprintf(stderr, "Complement Finished:%d\n", i);
	return;
} /* END */
