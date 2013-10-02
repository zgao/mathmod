#include <stdio.h>
#include <stdlib.h>

#include "heap.h"

inline int left(int i) { return 2 * i + 1; }
inline int right(int i) { return 2 * i + 2; }
inline int parent(int i) { return (i - 1) / 2; }

heap* make_heap(int maxsz) {
    int i;
    heap *ret = (heap*) malloc(sizeof(heap));
    ret->size = 0;
    ret->maxsz = maxsz;
    ret->id = (int*) malloc(sizeof(int) * maxsz);
    ret->lookup = (int*) malloc(sizeof(int) * maxsz);

    for (i = 0; i < maxsz; i++) {
        ret->id[i] = i;
        ret->lookup[i] = i;
    }

    ret->val = (double*) malloc(sizeof(double) * maxsz);
    return ret;
}

void destroy_heap(heap *h) {
    free( h->id );
    free( h->val );
    free( h->lookup );
    free( h );
}

int valid_heap(heap *h) {
    int *id = h->id, *lookup = h->lookup;
    int i;
    for (i = 0; i < h->size; i++) {
        if (!(0 <= id[i] && id[i] < h->maxsz)) return 0;
        if (lookup[id[i]] != i) return 0;
    }
    return 1;
}

void heap_swap(heap *h, int i, int j) {
    double s;
    int *lookup = h->lookup, *id = h->id;
    double *val = h->val;

    s = val[i];
    val[i] = val[j];
    val[j] = s;

    lookup[id[i]] = j;
    lookup[id[j]] = i;

    int s2;
    s2 = id[i];
    id[i] = id[j];
    id[j] = s2;
}

void heap_up(heap *h, int i) {
    if (i > 0 && h->val[parent(i)] > h->val[i]) {
        heap_swap(h, i, parent(i));
        heap_up(h, parent(i));
    }
}

void heap_down(heap *h, int i) {
    int a = left(i), b = right(i);
    double *val = h->val;
    if (b < h->size) {
        if (val[b] < val[a] && val[b] < val[i]) {
            heap_swap(h, i, b);
            heap_down(h, b);
            return;
        }
    }
    if (a < h->size && val[a] < val[i]) {
        heap_swap(h, i, a);
        heap_down(h, a);
    }
}

int empty(heap *h) {
    return h->size == 0;
}

pair pop(heap *h) {
    pair p;
    p.first = h->val[0];
    p.second = h->id[0];
    heap_swap(h, 0, h->size - 1);
    --(h->size);
    heap_down(h, 0);
    return p;
}

void update(heap *h, int q, double newval) {
    h->val[h->lookup[q]] = newval;
    heap_up(h, h->lookup[q]);
}
