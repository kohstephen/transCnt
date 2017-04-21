#include "infcylinder.h"

InfCylinder::InfCylinder(float radius, string mat, Temp t_init): Geometry(mat,t_init){
	_radius = radius;
}

InfCylinder::~InfCylinder(){}

float InfCylinder::getRadius(){
	return _radius;
}