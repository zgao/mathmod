#ifndef GENETIC_ALGORITH_H_ADDED
#define GENETIC_ALGORITH_H_ADDED

#include "GalleryArrangement.h"

void mutate(arrangement *x);
arrangement*  combine(arrangement *dad, arrangement *mom);
int selectGA(int *accumFitness, int length);
arrangement** generate(arrangement** previous, int length, float mutationRate, int elitism);

#endif
