#include "infcylinder.h"

InfCylinder::InfCylinder(Dim radius, string mat, Kelvin t_init): Geometry(mat,t_init){
	_radius = radius;
}

InfCylinder::~InfCylinder(){}

float InfCylinder::radius(){
	return _radius;
}