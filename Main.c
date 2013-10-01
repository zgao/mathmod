#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
	srand(time(NULL));
	arrangement *a = createRandomArrangement();
	arrangement *b = createRandomArrangement();
	printArrangement(a);
	printArrangement(b);
	printArrangement(combine(a,b));
	freeArrangement(a);
	freeArrangement(b);

	return 0;
}
