#include "GalleryArrangement.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

wall*  randomWall(int length, double startx, double startyi, double angleOffset);
void appendTo(wall *datum, wallList *x);

bool isAllowable(wallList *x){
	return true;
}

void freeArrangement(arrangement *x){
	wallList *current = x -> walls;
	while (current != NULL) {
		wallList *temp = current;
		freeWall(current -> value);
		current = current -> next;
		free(temp);
	}
	return;
}

void freeWall(wall *x){
	wall *current = x;
	while (current != NULL) {
		wall *temp = current;
		current = current -> child;
		free(temp);
	}
	return;
}

void appendTo(wall *datum, wallList *x) {
	if ( datum == NULL) {
		return;
	}
	if (x -> value == NULL) {
		x -> value = datum;
		return;
	} else {
		wallList *last = x;
		while(last -> next != NULL) {
			last = last -> next;
		}
		wallList *add = malloc(sizeof(wallList)); 
		add -> value = datum;
		add -> next = NULL;
		last -> next = add;
		return;
	}
}

wall* randomWall(int length, double startx, double starty, double angleOffset) {  //angleOffset is the angle subtracted from out -> angle to get the angle with respect to the vertical
	if (length == 0) {
		return NULL;
	} else {
		wall* out = malloc(sizeof(wall));
		out -> x_pos = startx;
		out -> y_pos = starty;
		double theta = M_PI*(double)rand()/(double)RAND_MAX;
		out -> angle = theta;
		double thetaV = theta - angleOffset;
		double newx = startx + 5.0*sin(thetaV);
		double newy = starty + 5*cos(thetaV);
		if (newx < 0 || newx > 22 || newy < 0 || newy >  20) {
			return NULL;
		} else {
			out -> child = randomWall(length - 1, newx, newy, M_PI_2 - thetaV);
			return out;
		}
	}
}

arrangement* createRandomArrangement() {
	arrangement *out = malloc(sizeof(arrangement));
	out -> cameraPattern = (int)rand() % 4;
	int numberOfWalls = 1 + (rand() % 5);
	int i;
	wallList *add;
	wallList *t = add;
	wallList *previous = NULL;
	for(i=0;i<numberOfWalls;i++) {
		add = malloc(sizeof(wallList));
		add -> value = randomWall(1 + rand() % 4, 22.0* (double)rand()/(double)RAND_MAX,20.0*(double)rand()/(double)RAND_MAX , 0.0);
		add -> next = NULL;
		if (previous != NULL) {
			previous -> next = add;
		}
		previous = add;
	}
	out -> walls = t;
	return out;
}

int wallLength(wall *x) {
	int out = 0;
	wall *curr = x;
	while(curr != NULL) {
		out += 1;
		curr = curr -> child;
	}
	return out;
}

double** wallPoints(wall *x) {
	int length = wallLength(x);
	double **out = malloc(length*sizeof(double*));
	int i = 0;
	wall *current = x;
	while (current != NULL) {
		out[i] = malloc(2*sizeof(double));
		out[i][0] = current -> x_pos;
		out[i][1] = current -> y_pos;
		i++;
	}
	return out;
}

int numberOfWalls(arrangement *x) {
	int out  = 0;
	wallList *current = x -> walls;
	while (current != NULL) {
		current = current -> next;
		out += 1;
	}
	return out;
}

void free2dArray(double **x, int length) {
	int i = 0;
	for(i= 0; i< length; i++) {
		free(x[i]);
	}
	free (x);
	return;
}

bool intersect(wall *a, wall *b) {
	double **apts = wallPoints(a);
	double **bpts = wallPoints(b);
	int aLen = wallLength(a);
	int bLen = wallLength(b);
	int i = 0;
	for (i = 0; i< aLen; i++) {
		int j = 0;
		for(j= 0; j < bLen; j++) {
			if(hypot(apts[i][0] - bpts[j][0], apts[i][1] - bpts[j][1])<= 4.9999999) {
				free2dArray(apts, aLen);
				free2dArray(bpts, bLen);
				return true;
			}
		}

	}
	free2dArray(apts, aLen);
	free2dArray(bpts, bLen);
	return false;
}
