#include "infcylinder.h"

InfCylinder::InfCylinder(Dim radius, string mat, Temp t_init): Geometry(mat,t_init){
	_radius = radius;
}

InfCylinder::~InfCylinder(){}

float InfCylinder::radius(){
	return _radius;
}