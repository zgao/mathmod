#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include "Camera.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
	srand(time(NULL));
	arrangement *a = createRandomArrangement();
	int i;
	for (i = 0; i < 100; i++) {
		mutate(a);
	}
	printArrangement(a);

	return 0;
}
