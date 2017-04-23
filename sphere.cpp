#include "sphere.h"

Sphere::Sphere(Dim radius, string mat, Temp t_init) : Geometry(mat,t_init){
		_radius = radius;
}

Sphere::~Sphere(){}

Dim Sphere::radius(){
	return _radius;
}
