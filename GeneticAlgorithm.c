#define WALL_PROBABILITY_TO_MUTATE 0.2
#define POSITION_PROBABILITY_TO_MUTATE 0.2
#define ANGLE_PROBABILITY_TO_MUTATE 0.1
#define CAMERA_PROBABILITY_TO_MUTATE 0.05

#include "GeneticAlgorithm.h"
#include "GalleryArrangement.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	double x;
	double y;
} point;

double BoxMuller(); //returns standard normal random number
wall* changeWall(wall *x);
point centerOfMass(wall *x);
int quadrant(point a);
void shuffle(arrangement **array, int n);
void weightedSort(arrangement **array, double *weights, int left, int right);
int  wsSeparate(arrangement **array, double *weights, int left, int right, int pivot);
void swapArr(arrangement **array, int a, int b);
void swapDoub(double *array, int a, int b);
void printWall(wall *x);

void printWall(wall *x) {
	if (x == NULL) {
		return;
	}
	printf("(%f,%f) -> ", x -> x_pos, x -> y_pos);
	printWall(x -> child);
}

void printArrangement(arrangement *x) {
	printf("ARRANGEMENT:\n");
	printf("Camera pattern:%d\n", x -> cameraPattern);
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
			int j = i + rand() / (RAND_MAX / (n-1) + 1);
			t = array[i];
			array[i] = array[j];
			array[j] = t;
		}
	}
	return;
}

point centerOfMass(wall *x) {
	int len = wallLength(x) + 1;
	double sumx = 0;
	double sumy = 0;
	wall *curr = x;
	while(curr != NULL) {
		sumx += curr -> x_pos;
		sumy += curr -> y_pos;
		curr = curr -> child;
	}
	return (point) {sumx/ (double) len, sumy / (double) len};
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
	if (c <= WALL_PROBABILITY_TO_MUTATE / 2 ) {
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
	}
	return x;
}

void mutate(arrangement *x){
	if ((double)rand() / (double) RAND_MAX < CAMERA_PROBABILITY_TO_MUTATE) {
		x -> cameraPattern += (rand() % 3 - 1);
	}
	wallList *current = x -> walls;
	while (current != NULL) {
		wall* mut = current -> value;
		if((double)rand() / (double)RAND_MAX < POSITION_PROBABILITY_TO_MUTATE) {
			double xchange = 0.5*BoxMuller();
			double ychange = 0.5*BoxMuller();
			while( mut != NULL ) {
				mut -> x_pos += xchange;
				mut -> y_pos += ychange;
				mut = mut -> child;
			}
		}
		mut = current -> value;
		while(mut != NULL){
			double angleChange = M_PI/9*BoxMuller();
			if((double)rand() / (double)RAND_MAX < ANGLE_PROBABILITY_TO_MUTATE && 0 <= (mut -> angle) && M_PI >= (mut -> angle)) {
				double centerX = mut -> x_pos;
				double centerY = mut -> y_pos;
				mut -> angle += angleChange;
				wall *sub = mut -> child;
				double r, theta;
				while( sub != NULL) {
					r = hypot((sub -> x_pos) - centerX, (sub -> y_pos) - centerY);
					theta = M_PI/2 + atan(((sub -> y_pos) - centerY) / ((sub -> x_pos) - centerX));
					sub -> x_pos = centerX + r*sin(theta + angleChange);
					sub -> y_pos = centerY + r*cos(theta + angleChange);
					sub = sub -> child;
				}
			}
			mut = mut -> child;
		}
		
		current -> value = changeWall(current -> value);

		current = current -> next;
	}
}

int quadrant(point a) {
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
	if ((double)rand() / RAND_MAX < 0.5) {
		out -> cameraPattern = dad -> cameraPattern;
	} else {
		out -> cameraPattern = mom -> cameraPattern;
	}
	out -> walls = malloc(sizeof(wallList));
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

arrangement** generate(arrangement **previous, int length, float mutationRate, int elitism, double (*fitnessp)(arrangement*)) {
	double *fitnesses = malloc(length*sizeof(double));
	arrangement **out = malloc(length*sizeof(arrangement*));
	int i;
	for(i = 0; i < length; i++) {
		fitnesses[i] = (*fitnessp)(previous[i]);
	}
	weightedSort(previous, fitnesses, 0, length - 1);
	int j = 0;
	for(i = length - 1; i > 0 ; i --) {
		for (j = i + 1; j < length ; j++) {
			fitnesses [j] += fitnesses[i];
		}
	}
	int parentPop = 2*(length - elitism);
	arrangement **parents = stochasticUniversalSample(previous, fitnesses, length, parentPop);
	shuffle(parents, parentPop);
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
