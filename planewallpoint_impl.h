#ifndef PLANEWALLPOINT_IMPL_H
#define PLANEWALLPOINT_IMPL_H
#include "constant.h"

class PlaneWallPoint_impl{
private:
	Loc _rect_loc;
	Secs _secs;
	Temp _temp;
public:
	PlaneWallPoint_impl(Loc rect_loc, Secs secs);
	~PlaneWallPoint_impl();
	Loc rect_loc();
	Temp temp();
	void temp(Temp temp);
	Secs time();
};

#endif