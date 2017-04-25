#include "cylinder.h"

Cylinder::Cylinder(Dim radius, Dim length, string mat, Kelvin t_init):Geometry(mat,t_init){
	_radius = radius;
	_length = length;
}

Cylinder::Cylinder(Dim radius, Dim length, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
	_radius = radius;
	_length = length;
}


Cylinder::~Cylinder(){}

Dim Cylinder::radius(){
	return _radius;
}

Dim Cylinder::length(){
	return _length;
}
