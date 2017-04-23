#include "planewallpoint.h"

PlaneWallPoint::PlaneWallPoint(Loc rect_loc, Secs secs) : Point(secs){
	_rect_loc = rect_loc;
}

PlaneWallPoint::~PlaneWallPoint(){}

Loc PlaneWallPoint::rect_loc(){
	return _rect_loc;
}

void PlaneWallPoint::rect_loc(Loc rect_loc) {
	_rect_loc = rect_loc;
}