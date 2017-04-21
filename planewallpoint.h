#ifndef PLANEWALLPOINT_IMPL_H
#define PLANEWALLPOINT_IMPL_H
#include "constant.h"

class PlaneWallPoint{
	Loc _rect_loc;
	Temp _temp;
	Secs _secs;
public:
	PlaneWallPoint(Loc rect_loc, Secs secs);
	~PlaneWallPoint();
	Loc rect_loc();
	Temp temp();
	void temp(Temp temp);
	Secs time();
};

#endif