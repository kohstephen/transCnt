#ifndef PLANEWALLPOINT_H
#define PLANEWALLPOINT_H
#include "point.h"

/**
 * PlaneWallPoint is a subclass of Point.
 * It is used to represent a point inside planewall.
 */
class PlaneWallPoint : public Point{
	Loc _rect_loc;

public:
	PlaneWallPoint(Loc rect_loc, Secs secs);
	~PlaneWallPoint();
	Loc rect_loc();
	void rect_loc(Loc rect_loc);
};

#endif