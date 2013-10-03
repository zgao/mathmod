#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include "Camera.h"
#include "geo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

inline double m2(double a, double b) {
    return (a > b) ? b : a;
}

inline double min(double a, double b, double c, double d) {
    return m2(m2(a, b), m2(c, d));
}

graph_wrapper atog(arrangement *a) {
    point *c = corners(a);
    point *p = paintings(a, c);
    return graph_of_arrangement(a, c, p);
}

int main(void) {
    srand(time(NULL));
    int i, j;
    //for (i = 0; i < 8; i++) {
        //printf("arrangement\n");

        arrangement *a = createRandomArrangement();
        printArrangement(a);
        graph_wrapper gw = atog(a);

        for (i = 0; i < 50; i++) {
            double *routes = shortest_paths(gw.graph, i, gw.size);
            printf("%d: %lf\n", i, min(routes[gw.size - 1], routes[gw.size - 2], routes[gw.size - 3], routes[gw.size - 4]));
        }


    //}

    //node *graph = graph_of_arrangement(a);

    return 0;
}
