#ifndef CAMERA_H_ADDED
#define CAMERA_H_ADDED
#include "GalleryArrangement.h"
#include <stdbool.h>

typedef struct {
	double x_p;
	double y_p;
} pt;

double angleFromCamera(double x, double y, int cam);
bool seenByCamera(double x, double y, int cam, arrangement *a);
pt wallIntersectSightLine(wallSequence *w, int cam, double x, double y);
double timeBetweenSight(bool sb1, bool sb2, double x, double y);


#endif
