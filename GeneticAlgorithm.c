#define WALL_PROBABILITY_TO_MUTATE 0.2
#define POSITION_PROBABILITY_TO_MUTATE 0.2
#define ANGLE_PROBABILITY_TO_MUTATE 0.1
#define CAMERA_PROBABILITY_TO_MUTATE 0.05

#include "GeneticAlgorithm.h"
#include "GalleryArrangement.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double BoxMuller(); //returns standard normal random number
wall* changeWall(wall *x);

double BoxMuller() {
	double u = (double)rand() / (double)RAND_MAX;
	double v = (double)rand() / (double)RAND_MAX;
	return sqrt(-2.0*log(u))*cos(2*M_PI*v);
}

wall* changeWall(wall *x) {
	double c = (double)rand() / (double)RAND_MAX;
	if (c <= WALL_PROBABILITY_TO_MUTATE / 2 ) {
		wall *out = malloc(sizeof(wall));
		double theta = (double)rand() / (double)RAND_MAX * M_PI;
		out -> child = x;
		out -> angle = theta;
		out -> x_pos = x -> x_pos - 5.0*sin(theta);
		out -> y_pos = x -> y_pos - 5.0*cos(theta);
		return out;
	} else if (c <= WALL_PROBABILITY_TO_MUTATE) {
		double theta = (double)rand() / (double)RAND_MAX * M_PI;
		wall *new = malloc(sizeof(wall));
		wall *last = x;
		while ( last -> child != NULL) {
			last = last -> child;
		}
		new -> child = NULL;
		last -> child = new;
		new -> angle = theta;
		new -> x_pos = last -> x_pos + 5.0*(sin(last -> angle));
		new -> y_pos = last -> y_pos + 5.0*(cos(last -> angle));
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
