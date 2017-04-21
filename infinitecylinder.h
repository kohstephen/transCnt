#ifndef INFINITECYLINDER_H
#define INFINITECYLINDER_H
#include "geometry.h"

class InfiniteCylinder : public Geometry{
	float _radius;
public:
	InfiniteCylinder(float radius, string mat, Temp t_init);
	~InfiniteCylinder();
	float getRadius();
};

#endif