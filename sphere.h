#ifndef SPHERE_H
#define SPHERE_H
#include "geometry.h"

class Sphere : public Geometry{
	float _radius;

public:
	Sphere(float radius, string mat, Temp t_init);
	~Sphere();
	float radius();
};

#endif