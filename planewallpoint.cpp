#include "planewallpoint.h"

PlaneWallPoint::PlaneWallPoint(Loc rect_loc, Secs secs){
	_rect_loc = rect_loc;
	_secs = secs;
	_temp = 0.0f;
}

PlaneWallPoint::~PlaneWallPoint(){}

Loc PlaneWallPoint::rect_loc(){
	return _rect_loc;
}

void PlaneWallPoint::temp(Temp temp) {
	_temp = temp;
}

Temp PlaneWallPoint::temp() {
	return _temp;
}

Secs PlaneWallPoint::time() {
	return _secs;
}