#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include <stdio.h>

int main(void) {
	arrangement *a = createRandomArrangement();
	mutate(a);
	freeArrangement(a);
	return 0;
}
