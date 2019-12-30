#include <stdlib.h>

void nulls(char line[], int n) {
	int i;
	for (i = 0; i < n + 1; i++)
		line[i] = '\00';
}
