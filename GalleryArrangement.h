#ifndef ARRANGEMENT_H_ADDED
#define ARRANGEMENT_H_ADDED

#include <deque>

using namespace std;

typedef struct wall {
    double x_pos; //define x position of the initial point of the wall
    double y_pos;
    double angle; //angle from 0 to pi of the wall, with 0 directly up
} wall;

typedef deque<wall*> wallSequence;
typedef deque<wallSequence> arrangement;

wallSequence* randomWall(int length, double startx, double starty, double angleOffset, deque<wall> *existing);
arrangement* createRandomArrangement();

#endif
