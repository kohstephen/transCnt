#include "cylinder.h"

Cylinder::Cylinder(Dim radius, Dim length, string mat, Temp t_init):Geometry(mat,t_init){
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
