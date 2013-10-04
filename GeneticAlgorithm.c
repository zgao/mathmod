#define WALL_PROBABILITY_TO_MUTATE 0.2
#define POSITION_PROBABILITY_TO_MUTATE 0.2
#define ANGLE_PROBABILITY_TO_MUTATE 0.1
#define RUNNING_MEAN 2.94
#define RUNNING_SD 0.56
#define NUMBER_OF_SEGMENTS 10

#include "GeneticAlgorithm.h"
#include "GalleryArrangement.h"
#include "Camera.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//#include <omp.h>
#include "geo.h"
#include "graph.h"

typedef struct {
	double x;
	double y;
} Point;

inline double m2(double a, double b) {
	return (a > b) ? b : a;
}

inline double min(double a, double b, double c, double d) {
	return m2(m2(a,b), m2(c,d));
}

double BoxMuller(); //returns standard normal random number
wall* changeWall(wall *x);
Point centerOfMass(wall *x);
int quadrant(Point a);
void shuffle(arrangement **array, int n);
void weightedSort(arrangement **array, double *weights, int left, int right);
int  wsSeparate(arrangement **array, double *weights, int left, int right, int pivot);
void quickSort(double *array, int left, int right);
int qsSeparate(double *array, int left, int right, int pivot);
void swapArr(arrangement **array, int a, int b);
void swapDoub(double *array, int a, int b);
void printWall(wall *x);
int f(double *securities, double speed);
double normalDF(double x);

double normalDF(double x) {
	return 1.0 / (RUNNING_SD * sqrt(2*M_PI)) * exp((-1.0) *  (x - RUNNING_MEAN)*(x - RUNNING_MEAN) / (2 * (RUNNING_SD) * (RUNNING_SD) ) );
}

int f(double *securities, double speed) {
	int j = 0;
	int i = 0;
	for(i = 0; i < 50; i ++) {
		if(securities[i] < speed) {
			j++;
		}
	}
	return j;
}

int compare(const void *a, const void *b) {
    return (*((double*)a) > *((double*)b));
}

double fitness(arrangement *a) {
        if (a->walls == NULL) return 0.0;
        wallList *wl;
        for (wl = a->walls; wl != NULL; wl = wl->next) {
            if (wl->value == NULL) return 0.0;
        }
        printf("%p\n", a->walls);
	point *corns = corners(a);
	point *paint = paintings(a, corns);
	if (paint == NULL) {
		//puts("PAINT IS NULL");
		return 0.0;
	}
	graph_wrapper G = graph_of_arrangement(a, paint, corns);
	double *securities = malloc(50*sizeof(double));
	int i;
	point *current = paint;
	for(i = 0; i < 50; i ++) {
		double *routes = shortest_paths(G.graph, i, G.size);
		double dist = min(routes[G.size - 1], routes[G.size - 2], routes[G.size - 3], routes[G.size - 4]);
		securities[i] = dist / (timeBetweenSight(seenByCamera(current -> x, current -> y, 1, a), seenByCamera(current -> x, current -> y, 2, a) , current -> x, current -> y));
	}
        qsort(securities, 50, sizeof(double), compare);
	double out = 0;
	double top = 6*RUNNING_SD + RUNNING_MEAN;
	double segmentLength = top / (double)NUMBER_OF_SEGMENTS;
	out += 0.5 * f(securities, 0.001) * normalDF(0.001);
	out += 0.5 * f(securities, top) * normalDF(0.001);
	for (i = 1; i < NUMBER_OF_SEGMENTS - 1; i ++) {
		out += f(securities, i*segmentLength) * normalDF(i*segmentLength);
	}
	out *= segmentLength;
	return out;
}

void printWall(wall *x) {
	if (x == NULL) {
		return;
	}
	printf("(%f,%f) -> ", x -> x_pos, x -> y_pos);
	printWall(x -> child);
}

void printArrangement(arrangement *x) {
	printf("ARRANGEMENT:\n");
	wallList *it = x -> walls;
	while(it != NULL) {
		printWall(it -> value);
		printf("\n");
		it = it -> next;
	}
}

void swapArr(arrangement **array, int a, int b) {
	arrangement *t = array[a];
	array[a] = array[b];
	array[b] = t;
}

void swapDoub(double *array, int a, int b) {
	double t = array[a];
	array[a] = array[b];
	array[b] = t;
}


void quickSort(double *array, int left, int right) {
	if ( left < right) {
		int piv = left+right / 2;
		int newPiv = qsSeparate(array, left, right, piv);
		quickSort(array, left, newPiv-1);
		quickSort(array, newPiv + 1, right);
	}	
}

int qsSeparate(double *array, int left, int right, int pivot) {
	double pivValue = array[pivot];
	swapDoub(array, pivot, right);
	int tempIndex = left;
	int i;
	for(i = left; i < right; i++) {
		if (array[i] < pivValue) {
			swapDoub(array, i, tempIndex);
			tempIndex ++;
		}
	}
	swapDoub(array, tempIndex, right );
}

int wsSeparate(arrangement **array, double *weights, int left, int right, int pivot) {
	double pivValue = weights[pivot];
	swapArr(array, pivot, right);
	swapDoub(weights, pivot, right);
	int tempIndex = left;
	int i;
	for(i = left; i < right; i++) {
		if (weights[i] > pivValue) {
			swapArr(array, i, tempIndex);
			swapDoub(weights, i, tempIndex);
			tempIndex ++;
		}
	}
	swapArr(array,tempIndex,right);
	swapDoub(weights, tempIndex, right );
}


void weightedSort(arrangement **array, double *weights, int left, int right) {
	if ( left < right) {
		int piv = left+right / 2;
		int newPiv = wsSeparate(array, weights, left, right, piv);
		weightedSort(array, weights, left, newPiv-1);
		weightedSort(array, weights, newPiv + 1, right);
	}	
}

void shuffle(arrangement **array, int n) {
	if ( n > 1 ) {
		int i;
		arrangement *t;
		for (i = 0; i < n-1; i++) {
			int j = rand() % (n-1);
			t = array[i];
			array[i] = array[j];
			array[j] = t;
		}
	}
	return;
}

Point centerOfMass(wall *x) {
	int len = wallLength(x) + 1;
	double sumx = 0;
	double sumy = 0;
	wall *curr = x;
	while(curr != NULL) {
		sumx += curr -> x_pos;
		sumy += curr -> y_pos;
		curr = curr -> child;
	}
	return (Point) {sumx/ (double) len, sumy / (double) len};
}

double BoxMuller() {
	double u = (double)rand() / (double)RAND_MAX;
	double v = (double)rand() / (double)RAND_MAX;
	return sqrt(-2.0*log(u))*cos(2*M_PI*v);
}

wall* changeWall(wall *x) {
	if(x == NULL) {
		return x;
	}
	double c = (double)rand() / (double)RAND_MAX;
	if (c <= WALL_PROBABILITY_TO_MUTATE / 2.0 ) {
		double theta = (double)rand() / (double)RAND_MAX * M_PI;
		double newX = x -> x_pos - 5.0*sin(M_PI_2+theta);
		double newY = x -> y_pos - 5.0*cos(M_PI_2+theta);
		if(fabs(newX - 11.0) > 11.0 || fabs(newY - 10.0) > 10.0) {
			return x;
		}
		wall *out = malloc(sizeof(wall));
		out -> child = x;
		out -> angle = theta;
		out -> x_pos = newX;
		out -> y_pos = newY;
		x -> angle -= (M_PI_2 - theta);
		return out;
	} else if (c <= WALL_PROBABILITY_TO_MUTATE) {
		double theta = (double)rand() / (double)RAND_MAX * M_PI;
		wall *last = x;
		double vAngle = M_PI_2;
		while ( last -> child != NULL) {
			vAngle += (last -> angle) - M_PI_2;
			last = last -> child;
		}
		double newX = last -> x_pos + 5.0*(sin(vAngle));
		double newY = last -> y_pos + 5.0*(cos(vAngle));
		if (fabs(newX - 11.0) > 11.0 || fabs(newY - 10.0) > 10.0) {
			return x;
		}
		wall *new = malloc(sizeof(wall));
		new -> child = NULL;
		last -> child = new;
		new -> angle = theta;
		new -> x_pos = newX;
		new -> y_pos = newY;
		return x;
	} else if (c <= WALL_PROBABILITY_TO_MUTATE * 2.0 && wallLength(x) > 2) {
		wall *last = x;
		wall *llast = last -> child;
		while ( llast -> child != NULL) {
			last = last -> child;
			llast = llast -> child;
		}
		freeWall(llast);
		last -> child = NULL;
	}
	return x;
}

void mutate(arrangement *x){
	wallList *current = x -> walls;
	while (current != NULL) {
		wall* mut = current -> value;
		if((double)rand() / (double)RAND_MAX < POSITION_PROBABILITY_TO_MUTATE) {
			double xchange = 0.5*BoxMuller();
			double ychange = 0.5*BoxMuller();
			bool allowMutate = true;
			
			while(mut != NULL && allowMutate){
				allowMutate = allowMutate && (fabs(mut -> x_pos + xchange - 11) < 11) && (fabs(mut -> y_pos + ychange - 10) < 10);
				mut = mut -> child;
			}
			mut = current -> value;
			if (allowMutate) {
				while( mut != NULL ) {
					mut -> x_pos += xchange;
					mut -> y_pos += ychange;
					mut = mut -> child;
				}
			}
		}
		mut = current -> value;
		while(mut != NULL){
			double angleChange = M_PI/9*BoxMuller();
			if((double)rand() / (double)RAND_MAX < ANGLE_PROBABILITY_TO_MUTATE && 0 <= (mut -> angle) && M_PI >= (mut -> angle)) {
				double centerX = mut -> x_pos;
				double centerY = mut -> y_pos;
				double r, theta;
				wall *sub = mut -> child;
				bool allowMutate2 = true;
				bool a, b;
				while(sub != NULL && allowMutate2) {
					r = hypot((sub -> x_pos) - centerX, (sub -> y_pos) - centerY);
					theta = M_PI/2 + atan(((sub -> y_pos) - centerY) / ((sub -> x_pos) - centerX));
					a = fabs(centerX + r*sin(theta + angleChange) - 11) < 11;
					b = fabs(centerY + r*cos(theta + angleChange) - 10) < 10;
					allowMutate2 = allowMutate2 && a && b;
					sub = sub -> child;
				}
				sub = mut -> child;
				if (allowMutate2) {
					mut -> angle += angleChange;
					while( sub != NULL) {
						r = hypot((sub -> x_pos) - centerX, (sub -> y_pos) - centerY);
						theta = M_PI/2 + atan(((sub -> y_pos) - centerY) / ((sub -> x_pos) - centerX));
						sub -> x_pos = centerX + r*sin(theta + angleChange);
						sub -> y_pos = centerY + r*cos(theta + angleChange);
						sub = sub -> child;
					}
				}
			}
			mut = mut -> child;
		}
		
		current -> value = changeWall(current -> value);

		current = current -> next;
	}
}

int quadrant(Point a) {
	if ( (a.x <= 11) && (a.y <= 10) ) {
		return 0;
	} else if ( (a.x > 11) && (a.y <= 10) ){
		return 1;
	} else if ((a.x <= 11) && (a.y > 10) ) {
		return 2;
	} else {
		return 3;
	}
}

arrangement* combine(arrangement *dad, arrangement *mom) {
	arrangement *out = malloc(sizeof(arrangement));
	wallList *W = malloc(sizeof(wallList));
	out -> walls = W;
	(out -> walls) -> value = NULL;
	(out -> walls) -> next = NULL;
	wallList *quadrants[4]; 
	int i;
	for ( i = 0; i < 4; i ++) {
		quadrants[i] = malloc(2*sizeof(wallList));
		quadrants[i][0].value =  NULL;
		quadrants[i][0].next = NULL;
		quadrants[i][1].value =  NULL;
		quadrants[i][1].next = NULL;
	}
	wallList *elemD = dad -> walls;
	wallList *elemM = mom -> walls;
	int q;
	while(elemD != NULL) {
		q = quadrant(centerOfMass(elemD -> value));
		appendTo(elemD -> value, &quadrants[q][0]);
		elemD = elemD -> next;
	}
	while(elemM != NULL) {
		q = quadrant(centerOfMass(elemM -> value));
		appendTo(elemM -> value, &quadrants[q][1]);
		elemM = elemM -> next;
	}
	for (i=0;i<4;i++) {  //need to add quadranty things to the combined one
		wallList *iterD = &quadrants[i][0];
		wallList *iterM = &quadrants[i][1];
		while(iterD != NULL && iterM != NULL) {
			if ((double)rand()/(double)RAND_MAX < 0.5) {
				appendTo(iterD -> value, out -> walls);
			} else {
				appendTo(iterM -> value, out -> walls);
			}
			iterD = iterD -> next;
			iterM = iterM -> next;
		}
		while(iterD != NULL) {
			if ((double)rand()/(double)RAND_MAX < 0.5) {
				appendTo(iterD -> value, out -> walls);
				iterD = iterD -> next;
			}
		}
		while(iterM != NULL) {
			if ((double)rand()/(double)RAND_MAX < 0.5) {
				appendTo(iterM -> value, out -> walls);
				iterM = iterM -> next;
			}
		}
	}
	return out;
}	

arrangement** stochasticUniversalSample(arrangement **population, double *accumFitness, int length, int keep) {
	double fitness = accumFitness[length-1];
	double dist = fitness/(double)keep;
	double start = dist * ((double)rand() / (double) RAND_MAX);
	double *points = malloc(keep*sizeof(double));
	int i;
	for(i=0;i<keep;i++) {
		points[i] = start + i*dist;
	}
	return rouletteWheelSelection(population, accumFitness, points, length, keep);
}

arrangement** rouletteWheelSelection(arrangement **population, double* accumFitness, double *points, int popLen, int ptsLen) {
	arrangement **out = malloc(ptsLen*sizeof(arrangement *));
	int i = 0;
	int j = 0;
	for(i=0;i<ptsLen;i++) {
		while(accumFitness[j] < points[i]) {
			j++;
		}
		out[i] = population[j];
	}
	free(points);
	return out;
}

int compareFitness(const void *a, const void *b) {
    pear
        aa = *(pear*)a,
        bb = *(pear*)b;
    return aa.fitness > bb.fitness;
}
arrangement** generate(arrangement **previous, int length, float mutationRate, int elitism) {
	double *fitnesses = malloc(length*sizeof(double));
	arrangement **out = malloc(length*sizeof(arrangement*));
	int i;
        int nonzero = 0;
	//#pragma omp parallel for
	for(i = 0; i < length; i++) {
            if (previous[i] == NULL) {
                printf("NULL PREV %d\n", i);
            }
            if (fitness(previous[i]) > 0.0) nonzero++;
		fitnesses[i] = 50.0 - fitness(previous[i]);
	}

        printf("nonzero: %d\n", nonzero);

        pear *af =
            malloc(length*sizeof(pear));
        for (i = 0; i < length; i++) {
            af[i].a = previous[i];
            af[i].fitness = fitnesses[i];
        }
        qsort(af, length, sizeof(pear), compareFitness);
        for (i = 0; i < length; i++)
            previous[i] = af[i].a;

	int j = 0;
	for(i = length - 1; i > 0 ; i --) {
		for (j = i + 1; j < length ; j++) {
			fitnesses [j] += fitnesses[i];
		}
	}
	int parentPop = 2*(length - elitism);
	arrangement **parents = stochasticUniversalSample(previous, fitnesses, length, parentPop);
	shuffle(parents, parentPop - 1);
	for(i = 0; i < elitism; i++) {
		out[i] = previous[i];
	}
	int pindex;
	for(i = elitism; i < length; i++ ) {
		pindex = 2*(i-elitism);
		out[i] = combine(parents[pindex], parents[pindex+1]);
		if((double)rand()/(double)RAND_MAX < mutationRate) {
			mutate(out[i]);
		}
	}
	free(fitnesses);
	for(i = 0; i< parentPop; i++) {
		freeArrangement(parents[i]);
	}
	free(parents);
	for(i = elitism; i< length; i++) {
		freeArrangement(previous[i]);
	}
	free(previous);
	return out;
}	
