#ifndef CYLINDER_H
#define CYLINDER_H
#include "geometry.h"

class Cylinder : public Geometry{
	float _radius;
	float _L;
public:
	Cylinder(float radius, float L, string mat, Temp t_init);
	~Cylinder();
	float getRadius();
	float getL();
};

#endif