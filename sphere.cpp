#include "sphere.h"

Sphere::Sphere(Dim radius, string mat, Kelvin t_init) : Geometry(mat,t_init){
		_radius = radius;
}

Sphere::~Sphere(){}

Dim Sphere::radius(){
	return _radius;
}
