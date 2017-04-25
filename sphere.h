#ifndef SPHERE_H
#define SPHERE_H
#include "geometry.h"

class Sphere : public Geometry{
	Dim _radius;

public:
	Sphere(Dim radius, string mat, Kelvin t_init);
	Sphere(Dim radius, float k, float c, float p, Kelvin t_init);
	~Sphere();
	Dim radius();
};

#endif