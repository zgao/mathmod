#ifndef GEO_H_ADDED
#define GEO_H_ADDED

#include <vector>

#include "graph.h"
#include "GalleryArrangement.h"

using namespace std;

typedef struct point {
    double x, y;
} point;

vector<point>* corners(arrangement *a);
vector<point>* paintings(arrangement *a, vector<point> *c);
vector<node>* graph_of_arrangement(arrangement *a, vector<point> *c, vector<point> *p);

#endif
