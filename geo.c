#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geo.h"

const double threshold = 2.0;

int ls_inter(double x1, double y1, double x2, double y2,
             double x3, double y3, double x4, double y4) {
    if (x1 == x3 && y1 == y3 || x2 == x3 && y2 == y3
            || x1 == x4 && y1 == y4 || x2 == x4 && y2 == y4) return false;
    x2 -= x1; x3 -= x1; x4 -= x1; y2 -= y1; y3 -= y1; y4 -= y1;
    double dist12 = sqrt(x2 * x2 + y2 * y2);
    double c = x2 / dist12;
    double s = y2 / dist12;
    double new_x = x3 * c + y3 * s;
    y3 = y3 * c - x3 * s; x3 = new_x;
    new_x = x4 * c + y4 * s;
    y4 = y4 * c - x4 * s; x4 = new_x;
    if (y3 < 0. && y4 < 0. || y3 >= 0. && y4 >= 0.) return false;
    double pos12 = x4 + (x3 - x4) * y4 / (y4 - y3);
    if (pos12 < 0. || pos12 > dist12) return false;
    return true;
}

int no_inter(arrangement *a, double x1, double y1, double x2, double y2) {
    wallList *wl;
    for (wl = a->walls; wl != NULL; wl = wl->next) {
        wall *first = wl->value;
        wall *second;
        for (second = first->child; second != NULL;
                first = first->child, second = second->child) {
            if (ls_inter(x1, y1, x2, y2,
                         first->x_pos, first->y_pos,
                         second->x_pos, second->y_pos))
                return 0;
        }
    }
    return 1;
}

int works(point *corners, double x, double y) {
    point *c;
    for (c = corners; c != NULL; c = c->next) {
        if (hypot(x - c->x, y - c->y) < threshold) {
            //printf("doesn't work: %lf %lf\n", x, y);
            return false;
        }
    }
    return true;
}

point* paintings(arrangement *a, point *c) {
    point *ret = (point*) malloc(sizeof(point) * 50);
    int paintings_placed = 0;

    wallList *wl = a->walls;
    for (wl = a->walls; wl != NULL; wl = wl->next) {
        wall *first = wl->value;
        wall *second = first->child;
        for (second = first->child; second != NULL;
                first = first->child, second = second->child) {
            double epsilon = 1e-2;
            double x_delta = first->y_pos - second->y_pos;
            double y_delta = second->x_pos - first->x_pos;
            double x_pos = (first->x_pos + second->x_pos) / 2.;
            double y_pos = (first->y_pos + second->y_pos) / 2.;
            ret[paintings_placed].x = x_pos + epsilon * x_delta;
            ret[paintings_placed++].y = y_pos + epsilon * y_delta;
            if (paintings_placed == 50) break;
            ret[paintings_placed].x = x_pos - epsilon * x_delta;
            ret[paintings_placed++].y = y_pos - epsilon * y_delta;
            if (paintings_placed == 50) break;
        }
    }

    double counters_x[4] = { 0.0, 1.0, 21.0, 22.0 };
    double counters_y[4] = { 19.0, 20.0, 0.0, 1.0 };
    double dx[4] = { 0.0, 2.0, -2.0, 0.0 };
    double dy[4] = { -2.0, 0.0, 0.0, 2.0 };
    int i = 0;

    int ctr = 0;

    while (paintings_placed < 50) {
        ctr++;

        if (ctr == 32) return NULL;

        i %= 4;
        if (works(c, counters_x[i], counters_y[i])) {
            ret[paintings_placed].x = counters_x[i];
            ret[paintings_placed++].y = counters_y[i];
            //printf("%lf %lf\n", counters_x[i], counters_y[i]);
        }
        counters_x[i] += dx[i];
        counters_y[i] += dy[i];
        i++;
    }

    for (i = 0; i < 50; i++) {
        ret[i].next = ret + (i + 1);
    }
    ret[49].next = NULL;

    return ret;
}

point* corners(arrangement *a) {
    int num_points = 0;
    point *ret = (point*) malloc(sizeof(point) * 100);
    int i;
    for (i = 0; i < 100; i++) {
        ret[i].x = ret[i].y = 0.0;
        ret[i].next = NULL;
    }
    wallList *wl;
    for (wl = a->walls; wl != NULL; wl = wl->next) {
        printf("a walls %p\n", a->walls);
        printf("value %p\n", wl->value);
        wall *first = wl->value;
        wall *second = first->child;
        ret[num_points].x = first->x_pos;
        ret[num_points++].y = first->y_pos;
        for (second = first->child; second != NULL;
                first = first->child, second = second->child) {
            ret[num_points].x = second->x_pos;
            ret[num_points++].y = second->y_pos;
        }
    }

    for (i = 0; i < num_points - 1; i++) ret[i].next = ret + i + 1;
    ret[num_points - 1].next = NULL;

    return ret;
}

graph_wrapper graph_of_arrangement(arrangement *a, point *c, point *p) {
    double doors[4][2];

    doors[0][0] = 0.0;
    doors[0][1] = 0.0;
    doors[1][0] = 2.0;
    doors[1][1] = 0.0;
    doors[2][0] = 20.0;
    doors[2][1] = 20.0;
    doors[3][0] = 22.0;
    doors[3][1] = 20.0;

    int n_corners = 0;
    point *corner;
    for (corner = c; corner != NULL; corner = corner->next) n_corners++;

    int n_vertices = 54 + n_corners;

    node* graph = (node*) malloc(sizeof(node) * n_vertices);
    int graph_size = 0;
    int i;
    for (i = 0; i < 50; i++)
        graph[graph_size] = make_new_node(graph_size++, 1, n_vertices);
    for (i = 0; i < n_corners; i++)
        graph[graph_size] = make_new_node(graph_size++, 0, n_vertices);
    for (i = 0; i < 4; i++)
        graph[graph_size] = make_new_node(graph_size++, 2, n_vertices);

    for (i = 0; i < 50; i++) {
        int j = 0;
        point *cj;
        for (cj = c; cj != NULL; cj = cj->next, j++) {
            if (no_inter(a,
                        p[i].x, p[i].y,
                        cj->x, cj->y)) {
                //printf("painting %d -> corner %d\n", i, j);
                make_new_edge(graph, i, j + 50,
                        hypot(
                            p[i].x - cj->x,
                            p[i].y - cj->y
                            ));
            }
        }
        for (j = 0; j < 4; j++) {
            if (no_inter(a,
                        p[i].x, p[i].y,
                        doors[j][0], doors[j][1])) {
                //printf("painting %d -> door %d\n", i, j);
                make_new_edge(graph, i, j + 50 + n_corners,
                        hypot(
                            p[i].x - doors[j][0],
                            p[i].y - doors[j][1]
                            ));
            }
        }
    }
    i = 0;
    point *ci;
    for (ci = c; ci != NULL; ci = ci->next, i++) {
        int j = 0;
        point *cj;
        for (cj = ci->next; cj != NULL; cj = cj->next, j++) {
            if (no_inter(a, ci->x, ci->y,
                        cj->x, cj->y)) {
                //printf("corner %d -> corner %d\n", i, j);
                make_new_edge(graph, i + 50, j + 50,
                        hypot(
                            ci->x - cj->x,
                            ci->y - cj->y
                            ));
            }
        }
        for (j = 0; j < 4; j++) {
            if (no_inter(a, ci->x, ci->y,
                        doors[j][0], doors[j][1])) {
                //printf("corner %d -> door %d\n", i, j);
                make_new_edge(graph, i + 50, j + 50 + n_corners,
                        hypot(
                            ci->x - doors[j][0],
                            ci->y - doors[j][1]
                            ));
            }
        }
    }

    graph_wrapper gw;
    gw.graph = graph;
    gw.size = n_vertices;

    return gw;
}
