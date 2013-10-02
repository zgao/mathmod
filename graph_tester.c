/*
#include <stdio.h>
#include <stdlib.h>

#include "graph.h"

FILE *fin, *fout;
int N, M;

int main() {
    fin = fopen("graph.in", "r");
    fout = fopen("graph.out", "w");
    fscanf(fin, "%d%d", &N, &M);
    node *graph = make_sample_graph(N);
    int i;
    for (i = 0; i < M; i++) {
        int tx, ty;
        double tc;
        fscanf(fin, "%d%d%lf", &tx, &ty, &tc);
        make_new_edge(graph, tx, ty, tc);
    }
    double *dists = shortest_paths(graph, 0, N);
    for (i = 0; i < N; i++)
        fprintf(fout, "%lf\n", dists[i]);

    fclose(fin);
    fclose(fout);
    free(graph);
    free(dists);

    return 0;
}*/
