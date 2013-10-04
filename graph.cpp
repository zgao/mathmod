#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "graph.h"

using namespace std;

#define infinity 1e9

node* make_new_node(int id, int flag, int sz) {
    node *n = (node*) malloc(sizeof(node));
    n->id = id;
    n->flag = flag;
    n->edges = new vector<node*>();
    n->weights = new vector<double>();
    return n;
}

void make_new_edge(vector<node*> *all_nodes, int i, int j, double weight) {
    vector<node*> all = *all_nodes;
    all[i]->edges->push_back(all[j]);
    all[j]->edges->push_back(all[i]);
    all[i]->weights->push_back(weight);
    all[j]->weights->push_back(weight);
}

double* shortest_paths(vector<node*> *all, int source, int sz) {
    vector<node*> all_nodes = *all;

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

    /*
    for (i = 0; i < sz; i++) {
        printf("%d:\n", i);
        int j;
        for (j = 0; j < all_nodes[i].num_edges; j++)
            printf("%d\n", all_nodes[i].edges[j].id);
    }*/

    while (!empty(pq)) {
        cpair p = pop(pq);
        int now = p.second;
        ret[now] = p.first;

        //printf("Got here: %d at dist %lf\n", p.second, p.first);

        node this_node = *(all_nodes[now]);
        seen[now] = 1;
        vector<node*> e = *(this_node.edges);
        for (i = 0; i < this_node.edges->size(); i++) {
            node next = *e[i];
            if (!seen[next.id]) {
                double dist = (*this_node.weights)[i];
                if (ret[now] + dist < ret[next.id]) {
                    ret[next.id] = ret[now] + dist;
                    //printf("updating: %d %lf\n", next.id, ret[next.id]);
                    update(pq, next.id, ret[next.id]);
                    assert(valid_heap(pq));
                }
            }
        }
    }

    destroy_heap(pq);

    return ret;
}

void destroy_graph(vector<node*> *gw) {
    vector<node*> unw = *gw;
    int i;
    for (i = 0; i < unw.size(); i++) {
        free(unw[i]->edges);
        free(unw[i]->weights);
    }
}
