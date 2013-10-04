#include "GalleryArrangement.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <deque>

using namespace std;

#define GLOBAL_MOD 8

inline wallSequence* randomWall(int length, double startx, double starty, double angleOffset);

wallSequence* randomWall(int length, double startx, double starty, double angleOffset, wallSequence *existing) {
    if (length == 0) {
        return NULL;
    } else {
        wall *out = (wall*) malloc(sizeof(wall));
        out->x_pos = startx;
        out->y_pos = starty;
        double theta = M_PI * (double) rand() / (double) RAND_MAX;
        out->angle = theta;
        double thetaV = theta - angleOffset;
        double newx = startx + 5.0 * sin(thetaV);
        double newy = starty + 5.0 * cos(thetaV);
        if (newx < 0 || newx > 22 || newy < 0 || newy > 20) {
            return NULL;
        } else {
            wallSequence *ret = randomWall(length - 1, newx, newy, M_PI_2 - thetaV);
            if (ret == NULL) ret = new wallSequence();
            ret->push_front(out);
            return ret;
        }
    }
}

inline wallSequence* randomWall(int length, double startx, double starty, double angleOffset) {
    return randomWall(length, startx, starty, angleOffset, new wallSequence());
}

arrangement* createRandomArrangement() {
    arrangement *out = new arrangement(); //out = add
    int numberOfWalls = 1 + (rand() % 5);
    for (int i = 0; i < numberOfWalls; i++) {
        wallSequence *addValue = NULL;
        while (!addValue || addValue->size() < 2) {
            if (addValue) delete addValue;
            addValue = randomWall(2 + rand() % 5, 22.0* (double)rand()/(double)RAND_MAX,20.0*(double)rand()/(double)RAND_MAX , 0.0);
        }
        out->push_back(*addValue);
    }
    return out;
}

