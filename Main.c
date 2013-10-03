#include "GalleryArrangement.h"
#include "GeneticAlgorithm.h"
#include "Camera.h"
#include "geo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


graph_wrapper atog(arrangement *a) {
    point *c = corners(a);
    point *p = paintings(a, c);
    return graph_of_arrangement(a, c, p);
}

int main(int argc, char **argv) {
    int n_generations = atoi(argv[1]);
    int generation_size = atoi(argv[2]);

    printf("Number of generations: %d\n", n_generations);
    printf("Generation size: %d\n", generation_size);

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
