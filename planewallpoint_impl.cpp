#include "planewallpoint_impl.h"

PlaneWallPoint_impl::PlaneWallPoint_impl(Loc rect_loc, Secs secs){
	_rect_loc = rect_loc;
	_secs = secs;
	_temp = 0.0f;
}

PlaneWallPoint_impl::~PlaneWallPoint_impl(){}

Loc PlaneWallPoint_impl::rect_loc(){
	return _rect_loc;
}

void PlaneWallPoint_impl::temp(Temp temp) {
	_temp = temp;
}

Temp PlaneWallPoint_impl::temp() {
	return _temp;
}

Secs PlaneWallPoint_impl::time() {
	return _secs;
}