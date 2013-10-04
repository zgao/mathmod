#define WALL_PROBABILITY_TO_MUTATE 0.2
#define POSITION_PROBABILITY_TO_MUTATE 0.2
#define ANGLE_PROBABILITY_TO_MUTATE 0.1
#define RUNNING_MEAN 2.94
#define RUNNING_SD 0.56
#define NUMBER_OF_SEGMENTS 10

#include "GeneticAlgorithm.h"
#include "GalleryArrangement.h"
#include "Camera.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>
#include "geo.h"
#include "graph.h"

using namespace std;

typedef struct {
    double x;
    double y;
} Point;

inline double m2(double a, double b) {
    return (a > b) ? b : a;
}

inline double min(double a, double b, double c, double d) {
    return m2(m2(a,b), m2(c,d));
}

double BoxMuller(); //returns standard normal random number
wallSequence* changeWall(wallSequence *x);
Point centerOfMass(wall *x);
int quadrant(Point a);
void shuffle(arrangement **array, int n);
void weightedSort(arrangement **array, double *weights, int left, int right);
int  wsSeparate(arrangement **array, double *weights, int left, int right, int pivot);
void quickSort(double *array, int left, int right);
int qsSeparate(double *array, int left, int right, int pivot);
void swapArr(arrangement **array, int a, int b);
void swapDoub(double *array, int a, int b);
void printWall(wall *x);
int f(double *securities, double speed);
double normalDF(double x);

double normalDF(double x) {
    return 1.0 / (RUNNING_SD * sqrt(2*M_PI)) * exp((-1.0) *  (x - RUNNING_MEAN)*(x - RUNNING_MEAN) / (2 * (RUNNING_SD) * (RUNNING_SD) ) );
}

int f(double *securities, double speed) {
    int j = 0;
    int i = 0;
    for(i = 0; i < 50; i ++) {
        if(securities[i] < speed) {
            j++;
        }
    }
    return j;
}

int compare(const void *a, const void *b) {
    return (*((double*)a) > *((double*)b));
}

bool intersect(wallSequence *a, wallSequence *b) {
    wallSequence as = *a, bs = *b;
    for (int i = 0; i < as.size() - 1; i++)
        for (int j = 0; j < bs.size() - 1; j++) {
            if (ls_inter(as[i]->x_pos, as[i]->y_pos, as[i + 1]->x_pos, as[i + 1]->y_pos,
                        bs[j]->x_pos, bs[j]->y_pos, bs[j + 1]->x_pos, bs[j + 1]->y_pos))
                return true;
        }
    return false;
}

bool selfIntersect(wallSequence *a) {
    wallSequence as = *a;
    for (int i = 0; i < as.size(); i++)
        for (int j = i + 1; j < as.size(); j++)
            if (hypot(as[i]->x_pos - as[j]->x_pos, as[i]->y_pos - as[j]->y_pos) < 3.54)
                return true;
    return false;
}

double fitness(arrangement *a) {
    arrangement pa = *a;
    if (!a) {
        return 0.0;
    }
    for (arrangement::iterator it = a->begin(); it != a->end(); it++) {
        wallSequence ws = *it;
        if (selfIntersect(&ws)) return 0.0;
        for (arrangement::iterator it2 = it + 1; it2 != a->end(); it2++) {
            if (intersect(&ws, &*it2)) {
                return 0.0;
            }
        }
    }
    vector<point> *c = corners(a);
    //printf("size of c is %d\n", (int) c->size());
    vector<point> *p = paintings(a, c);
    if (p == NULL) {
        return 0.0;
    }
    vector<point> corns = *c;
    
    for (int i = 0; i < corns.size(); i++)
        for (int j = i + 1; j < corns.size(); j++) {
            if (hypot(corns[i].x - corns[j].x, corns[i].y - corns[j].y) < 2.000)
                return 0.0;
        }
    vector<point> paint = *p;
    //printf("next size of c is %d\n", (int) c->size());
    vector<node*> *G = graph_of_arrangement(a, c, p);
    double *securities = (double*) malloc(50*sizeof(double));
    int i;
    for(i = 0; i < 50; i++) {
        double *routes = shortest_paths(G, i, G->size());
        double dist = min(routes[G->size() - 1],
                routes[G->size() - 2],
                routes[G->size() - 3],
                routes[G->size() - 4]);
        //printf("%lf\n", dist);
        point current = paint[i];
        securities[i] = dist / (timeBetweenSight(seenByCamera(current.x, current.y, 1, a), seenByCamera(current.x, current.y, 2, a) , current.x, current.y));
        free(routes);
    }
    qsort(securities, 50, sizeof(double), compare);
    double out = 0;
    double top = 6*RUNNING_SD + RUNNING_MEAN;
    double segmentLength = top / (double)NUMBER_OF_SEGMENTS;
    out += 0.5 * f(securities, 0.001) * normalDF(0.001);
    out += 0.5 * f(securities, top) * normalDF(0.001);
    for (i = 1; i < NUMBER_OF_SEGMENTS - 1; i ++) {
        out += f(securities, i*segmentLength) * normalDF(i*segmentLength);
    }
    out *= segmentLength;

    free(securities);
    destroy_graph(G);

    return out;
}

void printWall(wallSequence *x) {
    wallSequence xx = *x;
    for (int i = 0; i < xx.size(); i++) {
        printf("(%f,%f) -> ", xx[i]->x_pos, xx[i]->y_pos);
    }
    printf("\n");
}

void printArrangement(arrangement *x) {
    printf("ARRANGEMENT:\n");
    for (arrangement::iterator it = x->begin(); it != x->end(); it++) {
        printWall(&*it);
    }
}

void shuffle(arrangement **array, int n) {
    if ( n > 1 ) {
        int i;
        arrangement *t;
        for (i = 0; i < n-1; i++) {
            int j = rand() % (n-1);
            t = array[i];
            array[i] = array[j];
            array[j] = t;
        }
    }
    return;
}

Point centerOfMass(wallSequence *x) {
    double sumx = 0;
    double sumy = 0;
    for (wallSequence::iterator it = x->begin(); it != x->end(); it++) {
        wall *ref = *it;
        sumx += ref->x_pos;
        sumy += ref->y_pos;
    }
    double len = (double) (int) x->size();
    return (Point) {sumx/ (double) len, sumy / (double) len};
}

double BoxMuller() {
    double u = (double)rand() / (double)RAND_MAX;
    double v = (double)rand() / (double)RAND_MAX;
    return sqrt(-2.0*log(u))*cos(2*M_PI*v);
}

wallSequence* changeWall(wallSequence *xx) {
    if(xx == NULL) {
        return xx;
    }
    wallSequence x = *xx;
    double c = (double)rand() / (double)RAND_MAX;
    if (c <= WALL_PROBABILITY_TO_MUTATE / 2.0 ) {
        double theta = (double)rand() / (double)RAND_MAX * M_PI;
        double newX = x[0]->x_pos - 5.0*sin(M_PI_2+theta);
        double newY = x[0]->y_pos - 5.0*cos(M_PI_2+theta);
        if(fabs(newX - 11.0) > 11.0 || fabs(newY - 10.0) > 10.0) {
            return xx;
        }
        wall *out = (wall*) malloc(sizeof(wall));
        out->angle = theta;
        out->x_pos = newX;
        out->y_pos = newY;
        xx->push_front(out);
        x[0]->angle -= (M_PI_2 - theta);
        return xx;
    } else if (c <= WALL_PROBABILITY_TO_MUTATE) {
        double theta = (double)rand() / (double)RAND_MAX * M_PI;
        wall *last;
        double vAngle = M_PI_2;
        for (wallSequence::iterator it = x.begin(); it != x.end(); it++) {
            wall *w = *it;
            last = w;
            vAngle += (w->angle) - M_PI_2;
        }
        double newX = last->x_pos + 5.0*(sin(vAngle));
        double newY = last->y_pos + 5.0*(cos(vAngle));
        if (fabs(newX - 11.0) > 11.0 || fabs(newY - 10.0) > 10.0) {
            return xx;
        }
        wall *n = (wall*) malloc(sizeof(wall));
        n->angle = theta;
        n->x_pos = newX;
        n->y_pos = newY;
        xx->push_back(n);
        return xx;
    } else if (c <= WALL_PROBABILITY_TO_MUTATE * 2.0 && xx->size() > 2) {
        xx->pop_back();
    }
    return xx;
}

void mutate(arrangement *x){
    arrangement xx = *x;
    for (int i = 0; i < xx.size(); i++) {
        wallSequence mut = xx[i];

        if((double)rand() / (double)RAND_MAX < POSITION_PROBABILITY_TO_MUTATE) {
            double xchange = 0.5*BoxMuller();
            double ychange = 0.5*BoxMuller();

            bool allowMutate = false;

            for (wallSequence::iterator it = mut.begin();
                    it != mut.end();
                    it++) {
                wall *w = *it;
                if (fabs(w->x_pos + xchange - 11) >= 11 && fabs(w->y_pos + ychange - 10) >= 10) {
                    allowMutate = false;
                    break;
                }
            }

            if (allowMutate) {
                for (wallSequence::iterator it = mut.begin();
                        it != mut.end();
                        it++) {
                    wall *w = *it;
                    w->x_pos += xchange;
                    w->y_pos += ychange;
                }
            }
        }
        for (wallSequence::iterator it = mut.begin(); it != mut.end(); it++) {
            wall *w = *it;
            double angleChange = M_PI/9*BoxMuller();
            if((double)rand() / (double)RAND_MAX < ANGLE_PROBABILITY_TO_MUTATE && 0 <= (w->angle) && M_PI >= (w->angle)) {
                double centerX = w->x_pos;
                double centerY = w->y_pos;
                double r, theta;
                bool allowMutate2 = true;
                bool a, b;
                for (wallSequence::iterator sub = it + 1; sub != mut.end() && allowMutate2; sub++) {
                    wall *w2 = *sub;
                    r = hypot(w2->x_pos - centerX, w2->y_pos - centerY);
                    theta = M_PI/2 + atan((w->y_pos - centerY) / (w->x_pos - centerX));
                    a = fabs(centerX + r*sin(theta + angleChange) - 11) < 11;
                    b = fabs(centerY + r*cos(theta + angleChange) - 10) < 10;
                    allowMutate2 = allowMutate2 && a && b;
                }
                if (allowMutate2) {
                    for (wallSequence::iterator sub = it + 1; sub != mut.end(); sub++) {
                        wall *w2 = *sub;
                        w->angle += angleChange;
                        r = hypot(w->x_pos - centerX, w->y_pos - centerY);
                        theta = M_PI/2 + atan(w->y_pos - centerY) / (w->x_pos - centerX);
                        w->x_pos = centerX + r*sin(theta + angleChange);
                        w->y_pos = centerY + r*cos(theta + angleChange);
                    }
                }
            }
        }

        xx[i] = *changeWall(&(xx[i]));
    }
}

int quadrant(Point a) {
    if ( (a.x <= 11) && (a.y <= 10) ) {
        return 0;
    } else if ( (a.x > 11) && (a.y <= 10) ){
        return 1;
    } else if ((a.x <= 11) && (a.y > 10) ) {
        return 2;
    } else {
        return 3;
    }
}

arrangement* combine(arrangement *dad, arrangement *mom) {
    arrangement *out = new arrangement();
    arrangement *quadrants[4][2];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            quadrants[i][j] = new arrangement();
        }
    }
    int q;
    arrangement realDad = *dad, realMom = *mom;
    for (int i = 0; i < realDad.size(); i++) {
        wallSequence ws = realDad[i];
        q = quadrant(centerOfMass(&ws));
        quadrants[q][0]->push_back(ws);
    }
    for (int i = 0; i < realMom.size(); i++) {
        wallSequence ws = realMom[i];
        q = quadrant(centerOfMass(&ws));
        quadrants[q][1]->push_back(ws);
    }
    for (int i = 0; i < 4; i++) {  //need to add quadranty things to the combined one
        arrangement::iterator iterD = quadrants[i][0]->begin();
        arrangement::iterator iterM = quadrants[i][1]->begin();
        for (; iterD != quadrants[i][0]->end() && iterM != quadrants[i][1]->end(); iterD++, iterM++) {
            if ((double)rand()/(double)RAND_MAX < 0.5) {
                out->push_back(*iterD);
            } else {
                out->push_back(*iterM);
            }
        }
        for (; iterD != quadrants[i][0]->end(); iterD++) {
            if ((double)rand()/(double)RAND_MAX < 0.5) {
                out->push_back(*iterD);
            }
        }
        for (; iterM != quadrants[i][1]->end(); iterM++) {
            if ((double)rand()/(double)RAND_MAX < 0.5) {
                out->push_back(*iterM);
            }
        }
    }
    return out;
}

arrangement** stochasticUniversalSample(arrangement **population, double *accumFitness, int length, int keep) {
    double fitness = accumFitness[length-1];
    double dist = fitness/(double)keep;
    double start = dist * ((double)rand() / (double) RAND_MAX);
    double *points = (double*) malloc(keep*sizeof(double));
    int i;
    for(i=0;i<keep;i++) {
        points[i] = start + i*dist;
    }
    return rouletteWheelSelection(population, accumFitness, points, length, keep);
}

arrangement** rouletteWheelSelection(arrangement **population, double* accumFitness, double *points, int popLen, int ptsLen) {
    arrangement **out = (arrangement**) malloc(ptsLen*sizeof(arrangement *));
    int i = 0;
    int j = 0;
    for(i=0;i<ptsLen;i++) {
        while(accumFitness[j] < points[i]) {
            j++;
        }
        out[i] = population[j];
    }
    free(points);
    return out;
}

int compareFitness(const void *a, const void *b) {
    pear
        aa = *(pear*)a,
           bb = *(pear*)b;
    return aa.fitness > bb.fitness;
}
arrangement** generate(arrangement **previous, int length, float mutationRate, int elitism) {
    double *fitnesses = (double*) malloc(length*sizeof(double));
    arrangement **out = (arrangement**) malloc(length*sizeof(arrangement*));
    int i;
    int nonzero = 0;
    #pragma omp parallel for
    for(i = 0; i < length; i++) {
        if (previous[i] == NULL) {
            printf("NULL PREV %d\n", i);
        }
        if (fitness(previous[i]) > 0.0) nonzero++;
        fitnesses[i] = 50.0 - fitness(previous[i]);
    }

    printf("nonzero: %d\n", nonzero);

    pear *af =
        (pear*) malloc(length*sizeof(pear));
    for (i = 0; i < length; i++) {
        af[i].a = previous[i];
        af[i].fitness = fitnesses[i];
    }
    qsort(af, length, sizeof(pear), compareFitness);
    printf("Best: %lf\n", af[0].fitness);
    for (i = 0; i < length; i++)
        previous[i] = af[i].a;

    int j = 0;
    for(i = length - 1; i > 0 ; i --) {
        for (j = i + 1; j < length ; j++) {
            fitnesses [j] += fitnesses[i];
        }
    }
    int parentPop = 2*(length - elitism);
    arrangement **parents = stochasticUniversalSample(previous, fitnesses, length, parentPop);
    shuffle(parents, parentPop - 1);
    for(i = 0; i < elitism; i++) {
        out[i] = previous[i];
    }
    int pindex;
    for(i = elitism; i < length; i++ ) {
        pindex = 2*(i-elitism);
        out[i] = combine(parents[pindex], parents[pindex+1]);
        if((double)rand()/(double)RAND_MAX < mutationRate) {
            mutate(out[i]);
        }
    }
    free(fitnesses);
    free(parents);
    free(previous);
    return out;
}	
