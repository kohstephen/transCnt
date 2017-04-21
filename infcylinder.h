#ifndef INFCYLINDER_H
#define INFCYLINDER_H
#include "geometry.h"

class InfCylinder : public Geometry{
	float _radius;
public:
	InfCylinder(float radius, string mat, Temp t_init);
	~InfCylinder();
	float getRadius();
};

#endif