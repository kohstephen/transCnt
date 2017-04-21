#include "infinitecylinder.h"

InfiniteCylinder::InfiniteCylinder(float radius, string mat, Temp t_init): Geometry(mat,t_init){
	_radius = radius;
}

InfiniteCylinder::~InfiniteCylinder(){}

float InfiniteCylinder::getRadius(){
	return _radius;
}