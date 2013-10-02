#ifndef GEO_H_ADDED
#define GEO_H_ADDED

#include "graph.h"
#include "GalleryArrangement.h"

typedef struct {
    double x, y;
} point;

node* graph_of_arrangement(arrangement *a);

#endif
