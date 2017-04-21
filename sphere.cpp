#include "sphere.h"

Sphere::Sphere(float radius, string mat, Temp t_init) : Geometry(mat,t_init){
		_radius = radius;
}

Sphere::~Sphere(){}

float Sphere::radius(){
	return _radius;
}
