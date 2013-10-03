#ifndef GENETIC_ALGORITH_H_ADDED
#define GENETIC_ALGORITH_H_ADDED

#include "GalleryArrangement.h"

void mutate(arrangement *x);
arrangement*  combine(arrangement *dad, arrangement *mom);
arrangement** generate(arrangement** previous, int length, float mutationRate, int elitism);
arrangement** stochasticUniversalSample(arrangement **population, double *accumFitness, int length, int keep);
arrangement** rouletteWheelSelection(arrangement **population, double *accumFitness, double *points, int popLen, int ptsLen);
void printArrangement(arrangement *x);
double fitness(arrangement *a);

#endif
