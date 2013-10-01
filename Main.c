#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
	srand(time(NULL));
	arrangement *a = createRandomArrangement();
	printArrangement(a);
	int i;
	mutate(a);
	printArrangement(a);
	freeArrangement(a);

	//printWall(randomWall(1 + rand() % 4, 22.0* (double)rand()/(double)RAND_MAX, 20.0* (double)rand()/(double)RAND_MAX, 0.0));
	
	return 0;
}
