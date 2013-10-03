#ifndef GEO_H_ADDED
#define GEO_H_ADDED

#include "graph.h"
#include "GalleryArrangement.h"

typedef struct point {
    double x, y;
    struct point *next;
} point;

point* corners(arrangement *a);
point* paintings(arrangement *a, point *c);
graph_wrapper graph_of_arrangement(arrangement *a, point *c, point *p);

#endif
