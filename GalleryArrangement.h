#ifndef ARRANGEMENT_H_ADDED
#define ARRANGEMENT_H_ADDED
#include <stdbool.h>

typedef struct wall{
    double x_pos; //define x position of the initial point of the wall
    double y_pos;
    double angle; //angle from 0 to pi of the wall, with 0 directly up
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

wall*  randomWall(int length, double startx, double startyi, double angleOffset);
void appendTo(wall *datum, wallList *x);
bool isAllowable(wallList *x);
void freeArrangement(arrangement *x);
arrangement* createRandomArrangement();
void freeWall(wall *x);
int wallLength(wall *x);
int numberOfWalls(arrangement *x);
bool intersect(wall *x, wall *y);
void free2dArray(double **x, int length);

#endif
