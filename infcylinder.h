#ifndef INFCYLINDER_H
#define INFCYLINDER_H
#include "geometry.h"

class InfCylinder : public Geometry{
	Dim _radius;
public:
	InfCylinder(Dim radius, string mat, Kelvin t_init);
	~InfCylinder();
	Dim radius();
};

#endif