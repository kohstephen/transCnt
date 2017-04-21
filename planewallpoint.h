#ifndef PLANEWALLPOINT_H
#define PLANEWALLPOINT_H
#include "point.h"

class PlaneWallPoint : public Point{
	Loc _rect_loc;

public:
	PlaneWallPoint(Loc rect_loc, Secs secs);
	~PlaneWallPoint();
	Loc rect_loc();
};

#endif