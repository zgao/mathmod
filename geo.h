#ifndef GEO_H_ADDED
#define GEO_H_ADDED

#include <vector>

#include "graph.h"
#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"

using namespace std;

typedef struct point {
    double x, y;
} point;

int ls_inter(double x1, double y1, double x2, double y2,
        double x3, double y3, double x4, double y4);

vector<point>* corners(arrangement *a);
vector<point>* paintings(arrangement *a, vector<point> *c);
vector<node*>* graph_of_arrangement(arrangement *a, vector<point> *c, vector<point> *p);

#endif
