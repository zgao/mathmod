#include "Camera.h"
#include "GalleryArrangement.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double angleFromCamera(double x, double y, int cam) { //measured from horizontal in obvious way
	if (cam == 1) {
		return angleFromCamera(22 - x, 20 - y, 2);
	}
	return atan(y/(22-x));
}

bool seenByCamera(double x, double y, int cam, arrangement *a) {
	bool p = cam == 1;
	bool q = hypot(20 - y, x) < 4.47;
	bool r = hypot(22 - x, y) < 4.47;
	if ( (p && q) || (!p && r) ) {
		return false;
	}
	if (p) {
		wallList * curr = a -> walls;
		double closestX = 0.1;
		double closestY = 19.9;
		double dist = hypot(closestX, 20 -closestY);
		while (curr != NULL && curr -> value != NULL) {
			pt inter = wallIntersectSightLine(curr -> value, cam, x, y);
			if (inter.x_p > 0.0 && hypot(x-inter.x_p, y-inter.y_p) < hypot(x-closestX, y-closestY) ) {
				closestX = inter.x_p;
				closestY = inter.y_p;
			}
			curr = curr -> next;
		}
		return ( hypot(x, 20.0-y) >= 3.0*hypot(closestX, 20.0- closestY)) ;
	} else {
		wallList * curr = a -> walls;
		double closestX = 21.9;
		double closestY = 0.1;
		double dist = hypot(22.0 - closestX, closestY);
		while (curr != NULL && curr -> value != NULL) {
			pt inter = wallIntersectSightLine(curr -> value, cam, x, y);
			if (inter.x_p > 0.0 && hypot(x-inter.x_p, y-inter.y_p) < hypot(x-closestX, y-closestY) ) {
				closestX = inter.x_p;
				closestY = inter.y_p;
			}
			curr = curr -> next;
		}
		return ( hypot(22.0 - x, y) >= 3.0*hypot(22.0 - closestX, closestY)) ;
	}
}

pt wallIntersectSightLine(wall *w, int cam, double x, double y) {
	if (w == NULL) {
		pt out = {-1.0, -1.0};
		return out;
	}
	double Cx, Cy;
	if (cam == 1) {
		Cx = 0.0;
		Cy = 20.0;
	} else {
		Cx = 22.0;
		Cy = 0.0;
	}
	double slope = (y - Cy) / (x - Cx);
	double intercept = y + slope*x;
	wall *s = w;
	wall *t = w -> child;
	while(t != NULL) {
		double m = ( (s -> y_pos) - (t -> y_pos) ) / ( (s -> x_pos) - (t -> x_pos) ) ;
		double b = (s-> y_pos) + m*(s -> x_pos);
		double intY = (intercept - b*slope/m)/(1 - slope / m );
		double intX = (intY - b) / m;
		if (hypot(Cx - intX, Cy - intY) < hypot ( Cx - x,Cy - y) ) {
			pt out = {intX, intY};
			return out;
		}
		wall * temp = t -> child;
		s = t;
		t = temp;
	}
	pt out = {-1.0, -1.0};
	return out;
}

double timeBetweenSight(bool sb1, bool sb2, double x, double y) {
	bool closerTo1 = (y > x*20.0/22.0);
	if(sb1 && ! sb2) {
		double a = angleFromCamera(x, y, 1);
		return 2*(fabs(M_PI_4 - a) + M_PI_4 - M_PI / 6.0) / (M_PI/30.0);
	}
	if(sb2 && ! sb1) {
		double a = angleFromCamera(x, y, 2);
		return 2*(fabs(M_PI_4 - a) + M_PI_4 - M_PI / 6.0) / (M_PI/30.0);
	}
	double a1 = angleFromCamera(x, y, 1);
	double a2 = angleFromCamera(x,y,2);
	bool mt1 = fabs(a1 - M_PI_4) < M_PI / 12.0;
	bool mt2 = fabs(a2 - M_PI_4) < M_PI / 12.0;
	if(closerTo1) {
		if (  mt1 && mt2){
			return 2 * (M_PI_4 - fabs(a2 - M_PI_4) - M_PI / 6.0 ) / (M_PI / 30.0);
		}
		if (mt2 && !mt1) {
			return fmax( 2*(M_PI /12.0 - fabs(a2-M_PI_4)) / (M_PI / 30.0) ,  (fabs(a1-M_PI_4) - M_PI / 12.0 - (M_PI /12.0 - fabs(a2-M_PI_4))) / (M_PI / 30.0));
		}
		if (!mt2 && mt1) {
			return fmax( 2*(M_PI /12.0 - fabs(a1-M_PI_4)) / (M_PI / 30.0) ,  (fabs(a2-M_PI_4) - M_PI / 12.0 - (M_PI /12.0 - fabs(a1-M_PI_4))) / (M_PI / 30.0));
		}
		return 2 * a2 / (M_PI / 30.0);
	} else {
		return timeBetweenSight(true, true, -x + 22.0, -y + 20.0);
	}
}
