#ifndef GRAPH_H_ADDED
#define GRAPH_H_ADDED

#include "heap.h"

typedef struct node {
    int id;
    int flag;

    int num_edges;
    struct node *edges; //dest nodes
    double *weights; //respective weights
} node;

typedef struct {
    node *graph;
    int size;
} graph_wrapper;

node* make_sample_graph(int sz);
node make_new_node(int id, int flag, int sz);
void make_new_edge(node *all, int i, int j, double weight);
double* shortest_paths(node *all_nodes, int source, int sz);
void destroy_graph(graph_wrapper gw);

#endif
