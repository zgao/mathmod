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
int compareFitness2(const void *a, const void *b) {
    pear
        aa = *(pear*)a,
        bb = *(pear*)b;
    return aa.fitness > bb.fitness;
}

int main(int argc, char **argv) {
    int n_generations = atoi(argv[1]);
    int generation_size = atoi(argv[2]);

    printf("Number of generations: %d\n", n_generations);
    printf("Generation size: %d\n", generation_size);

    srand(time(NULL));

    printArrangement(createRandomArrangement());
	
    arrangement **old = malloc(generation_size*sizeof(arrangement*));
    int i;
    for(i = 0; i < generation_size; i++) {
        old[i] = createRandomArrangement();
    }
    for(i = 0; i < n_generations; i++) {
        arrangement **new = generate(old, generation_size, 0.05, generation_size / 100 );
        int j;
        for(j = 0; j < generation_size; j++) {
            //freeArrangement(old[j]);
        }
        //free(old);
        old = new;
    }
    double *fitnesses = malloc(generation_size*(sizeof(double)));
    for(i = 0; i < generation_size; i ++) {
        fitnesses[i] = 50 - fitness(old[i]);
    }

    pear *af = malloc(generation_size*sizeof(pear));
    for (i = 0; i < generation_size; i++) {
        af[i].a = old[i];
        af[i].fitness = fitnesses[i];
    }
    qsort(af, generation_size, sizeof(pear), compareFitness2);
    printArrangement(af[0].a);
    //int i, j;
    //for (i = 0; i < 8; i++) {
    //printf("arrangement\n");

    // arrangement *a = createRandomArrangement();
    // printArrangement(a);
    // graph_wrapper gw = atog(a);

    //for (i = 0; i < 50; i++) {
    //   double *routes = shortest_paths(gw.graph, i, gw.size);
    //    printf("%d: %lf\n", i, min(routes[gw.size - 1], routes[gw.size - 2], routes[gw.size - 3], routes[gw.size - 4]));
    //}


    //}

    //node *graph = graph_of_arrangement(a);

    return 0;
}
