#ifndef ARRANGEMENT_H_ADDED
#define ARRANGEMENT_H_ADDED
#include <stdbool.h>

typedef struct wall{
	int x_pos; //define x position of the initial point of the wall
	int y_pos;
	int angle; //angle from 0 to 180 of the child, with 0 directly up
        struct wall* child; 
} wall;

//we don't know how many walls we'll have, so we use a linked list of them
typedef struct wallList{
	wall* value;
	struct wallList* next;
} wallList;

typedef struct arrangement{
	wallList* walls;
	int cameraPattern; // we have a finite number of camera patterns we're testing
} arrangement;

bool isAllowable(wallList *x);
void freeArrangement(arrangement *x);
arrangement* createRandomArrangement();
void freeWall(wall *x);

#endif
