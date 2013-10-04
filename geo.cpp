#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "geo.h"

using namespace std;

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
    for (arrangement::iterator it = a->begin(); it != a->end(); it++) {
        wallSequence ws = *it;
        for (int i = 0; i < ws.size() - 1; i++) {
            if (ls_inter(x1, y1, x2, y2,
                        ws[i]->x_pos, ws[i]->y_pos, ws[i + 1]->x_pos, ws[i + 1]->y_pos))
                return 0;
        }
    }
    return 1;
}

int works(vector<point> *corners, double x, double y) {
    for (vector<point>::iterator it = corners->begin(); it != corners->end(); it++) {
        if (hypot(x - it->x, y - it->y) < threshold) return false;
    }
    return true;
}

vector<point>* paintings(arrangement *a, vector<point> *c) {
    vector<point> *ret = new vector<point>();
    int paintings_placed = 0;

    for (arrangement::iterator wp = a->begin(); wp != a->end(); wp++) {
        wallSequence ws = *wp;
        for (int i = 0; i < ws.size() - 1; i++) {
            wall *first = ws[i], *second = ws[i + 1];
            double epsilon = 1e-2;
            double x_delta = first->y_pos - second->y_pos;
            double y_delta = second->x_pos - first->x_pos;
            double x_pos = (first->x_pos + second->x_pos) / 2.;
            double y_pos = (first->y_pos + second->y_pos) / 2.;
            //printf("2 paintings placed %lf %lf\n", x_pos, y_pos);
            point new1, new2;
            new1.x = x_pos + epsilon * x_delta;
            new1.y = y_pos + epsilon * y_delta;
            ret->push_back(new1);
            if (ret->size() == 50) break;
            new2.x = x_pos - epsilon * x_delta;
            new2.y = y_pos - epsilon * y_delta;
            ret->push_back(new2);
            if (ret->size() == 50) break;
        }
    }

    double counters_x[4] = { 0.0, 1.0, 21.0, 22.0 };
    double counters_y[4] = { 19.0, 20.0, 0.0, 1.0 };
    double dx[4] = { 0.0, 2.0, -2.0, 0.0 };
    double dy[4] = { -2.0, 0.0, 0.0, 2.0 };
    int i = 0;

    int ctr = 0;

    while (ret->size() < 50) {
        ctr++;

        if (ctr == 36) {
            //printf("returning null now\n");
            return NULL;
        }

        i %= 4;
        if (works(c, counters_x[i], counters_y[i])) {
            point next;
            next.x = counters_x[i];
            next.y = counters_y[i];
            ret->push_back(next);
            //printf("%lf %lf\n", counters_x[i], counters_y[i]);
        }
        counters_x[i] += dx[i];
        counters_y[i] += dy[i];
        i++;
    }

    return ret;
}

vector<point>* corners(arrangement *a) {
    int num_points = 0;
    vector<point> *ret = new vector<point>();
    for (arrangement::iterator it = a->begin(); it != a->end(); it++) {
        wallSequence ws = *it;
        for (wallSequence::iterator wsit = ws.begin(); wsit != ws.end(); wsit++) {
            wall *w = *wsit;
            point p;
            p.x = w->x_pos;
            p.y = w->y_pos;
        }
    }
    return ret;
}

vector<node>* graph_of_arrangement(arrangement *a, vector<point> *cc, vector<point> *pp) {
    vector<point> p = *pp;
    vector<point> c = *cc;

    double doors[4][2];

    doors[0][0] = 0.0;
    doors[0][1] = 0.0;
    doors[1][0] = 2.0;
    doors[1][1] = 0.0;
    doors[2][0] = 20.0;
    doors[2][1] = 20.0;
    doors[3][0] = 22.0;
    doors[3][1] = 20.0;

    int n_corners = c.size();

    int n_vertices = 54 + n_corners;

    vector<node> *graph = new vector<node>();
    for (int i = 0; i < 50; i++)
        graph->push_back(make_new_node(graph->size(), 1, n_vertices));
    for (int i = 0; i < n_corners; i++)
        graph->push_back(make_new_node(graph->size(), 0, n_vertices));
    for (int i = 0; i < 4; i++)
        graph->push_back(make_new_node(graph->size(), 2, n_vertices));

    for (int i = 0; i < 50; i++) {
        int j = 0;
        for (vector<point>::iterator cj = c.begin(); cj != c.end(); cj++) {
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
    for (int i = 0; i < c.size(); i++) {
        vector<point>::iterator ci = c.begin() + i;
        for (int j = i + 1; j < c.size(); j++) {
            vector<point>::iterator cj = c.begin() + j;
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
        for (int j = 0; j < 4; j++) {
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

    return graph;
}
