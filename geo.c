#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geo.h"

const int n_paintings = 50;
const int n_corners = 150;
const int n_vertices = 200;
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

int works(point *corners, int nc, double x, double y) {
    int i;
    for (i = 0; i < nc; i++) {
        if (hypot(x - corners[i].x, y - corners[i].y) < threshold) {
            //printf("doesn't work: %lf %lf\n", x, y);
            return false;
        }
    }
    return true;
}

graph_wrapper graph_of_arrangement(arrangement *a) {
    double paintings[n_paintings][2]; //use epsilon normal
    double doors[4][2];

    doors[0][0] = 0.0;
    doors[0][1] = 0.0;
    doors[1][0] = 2.0;
    doors[1][1] = 0.0;
    doors[2][0] = 20.0;
    doors[2][1] = 20.0;
    doors[3][0] = 22.0;
    doors[3][1] = 20.0;

    int paintings_placed = 0, num_points = 0;
    point *corners = (point*) malloc(sizeof(point) * n_corners);
    //0,0 lower left hand corner, 20x22, 2 meter wide door
    wallList *wl;
    for (wl = a->walls; wl != NULL; wl = wl->next) {
        wall *first = wl->value;
        wall *second;
        corners[num_points].x = first->x_pos;
        corners[num_points++].y = first->y_pos;
        for (second = first->child; second != NULL;
                first = first->child, second = second->child) {
            corners[num_points].x = second->x_pos;
            corners[num_points++].y = second->y_pos;
            double epsilon = 1e-2;
            double x_delta = first->y_pos - second->y_pos;
            double y_delta = second->x_pos - first->x_pos;
            double x_pos = (first->x_pos + second->x_pos) / 2.;
            double y_pos = (first->y_pos + second->y_pos) / 2.;
            paintings[paintings_placed][0] = x_pos + epsilon * x_delta;
            paintings[paintings_placed++][1] = y_pos + epsilon * y_delta;
            paintings[paintings_placed][0] = x_pos - epsilon * x_delta;
            paintings[paintings_placed++][1] = y_pos - epsilon * y_delta;
            //printf("2 paintings: %lf %lf\n", x_pos, y_pos);
        }
    }
    double counters_x[4] = { 0.0, 1.0, 21.0, 22.0 };
    double counters_y[4] = { 19.0, 20.0, 0.0, 1.0 };
    double dx[4] = { 0.0, 2.0, -2.0, 0.0 };
    double dy[4] = { -2.0, 0.0, 0.0, 2.0 };
    int i = 0;
    while (paintings_placed < 50) {
        i %= 4;
        if (works(corners, num_points, counters_x[i], counters_y[i])) {
            paintings[paintings_placed][0] = counters_x[i];
            paintings[paintings_placed++][1] = counters_y[i];
            printf("%lf %lf\n", counters_x[i], counters_y[i]);
        }
        counters_x[i] += dx[i];
        counters_y[i] += dy[i];
        //check out of bounds
        i++;
    }
    //now that we have the painting positions, we can form the graph
    node* graph = (node*) malloc(sizeof(node) * n_vertices);
    int graph_size = 0;
    for (i = 0; i < n_paintings; i++)
        graph[graph_size] = make_new_node(graph_size++, 1, n_vertices);
    for (i = 0; i < num_points; i++)
        graph[graph_size] = make_new_node(graph_size++, 0, n_vertices);
    for (i = 0; i < 4; i++)
        graph[graph_size] = make_new_node(graph_size++, 2, n_vertices);
    for (i = 0; i < n_paintings; i++) {
        int j;
        for (j = 0; j < num_points; j++) {
            int flag = 0;
            if (no_inter(a,
                        paintings[i][0], paintings[i][1],
                        corners[j].x, corners[j].y)) {
                //printf("painting %d -> corner %d\n", i, j);
                make_new_edge(graph, i, j + n_paintings,
                        hypot(
                            paintings[i][0] - corners[j].x,
                            paintings[i][1] - corners[j].y
                            ));
            }
        }
        for (j = 0; j < 4; j++) {
            if (no_inter(a,
                        paintings[i][0], paintings[i][1],
                        doors[j][0], doors[j][1])) {
                //printf("painting %d -> door %d\n", i, j);
                make_new_edge(graph, i, j + n_paintings + num_points,
                        hypot(
                            paintings[i][0] - doors[j][0],
                            paintings[i][1] - doors[j][1]
                            ));
            }
        }
    }
    for (i = 0; i < num_points; i++) {
        int j;
        for (j = 0; j < num_points; j++) {
            if (i == j) continue;
            if (no_inter(a, corners[i].x, corners[i].y,
                        corners[j].x, corners[j].y)) {
                //printf("corner %d -> corner %d\n", i, j);
                make_new_edge(graph, i + n_paintings, j + n_paintings,
                        hypot(
                            corners[i].x - corners[j].x,
                            corners[i].y - corners[j].y
                            ));
            }
        }
        for (j = 0; j < 4; j++) {
            if (no_inter(a, corners[i].x, corners[i].y,
                        doors[j][0], doors[j][1])) {
                //printf("corner %d -> door %d\n", i, j);
                make_new_edge(graph, i + n_paintings, j + n_paintings + num_points,
                        hypot(
                            corners[i].x - doors[j][0],
                            corners[i].y - doors[j][1]
                            ));
            }
        }
    }

    graph_wrapper gw;
    gw.graph = graph;
    gw.size = n_paintings + num_points + 4;

    return gw;
}
