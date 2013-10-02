#ifndef HEAP_H_ADDED
#define HEAP_H_ADDED

#include "pair.h"

typedef struct {
    int size, maxsz;
    int *id;
    int *lookup;
    double *val;
} heap;

heap* make_heap     (int maxsz);
void destroy_heap   (heap *h);

int valid_heap      (heap *h);
void heap_swap      (heap *h, int i, int j);
void heap_up        (heap *h, int i);
void heap_down      (heap *h, int i);
int empty           (heap *h);
pair pop            (heap *h);
void update         (heap *h, int q, double newval);

#endif
