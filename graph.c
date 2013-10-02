#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "graph.h"

#define infinity 1e9

node* make_sample_graph(int sz) {
    node *ret = (node*) malloc(sizeof(node) * sz);
    int i;
    for (i = 0; i < sz; i++)
        ret[i] = make_new_node(i, 1, sz);
    return ret;
}

node make_new_node(int id, int flag, int sz) {
    node n;
    n.id = id;
    n.flag = flag;
    n.num_edges = 0;
    n.edges = (node*) malloc(sizeof(node) * sz);
    n.weights = (double*) malloc(sizeof(node) * sz);
    return n;
}

void make_new_edge(node *all, int i, int j, double weight) {
    int sz_i = all[i].num_edges++, sz_j = all[j].num_edges++;
    all[i].edges[sz_i] = all[j];
    all[j].edges[sz_j] = all[i];
    all[i].weights[sz_i] = weight;
    all[j].weights[sz_j] = weight;
}

double* shortest_paths(node *all_nodes, int source, int sz) {
    heap *pq = make_heap(sz);
    pq->size = sz;
    double *ret = (double*) malloc(sizeof(double) * sz);
    int *seen = (int*) malloc(sizeof(int) * sz);

    int i;
    for (i = 0; i < sz; i++) {
        ret[i] = infinity;
        pq->val[i] = infinity;
        seen[i] = 0;
    }

    update(pq, source, 0.0);
    ret[source] = 0.0;
    seen[source] = 1;

    while (!empty(pq)) {
        pair p = pop(pq);
        int now = p.second;
        ret[now] = p.first;
        node this_node = all_nodes[now];
        seen[now] = 1;
        for (i = 0; i < this_node.num_edges; i++) {
            node next = this_node.edges[i];
            if (!seen[next.id]) {
                double dist = this_node.weights[i];
                if (ret[now] + dist < ret[next.id]) {
                    ret[next.id] = ret[now] + dist;
                    update(pq, next.id, ret[next.id]);
                    assert(valid_heap(pq));
                }
            }
        }
    }

    destroy_heap(pq);

    return ret;
}
