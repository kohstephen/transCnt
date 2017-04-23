#ifndef INFCYLINDER_H
#define INFCYLINDER_H
#include "geometry.h"

class InfCylinder : public Geometry{
	Dim _radius;
public:
	InfCylinder(Dim radius, string mat, Temp t_init);
	~InfCylinder();
	Dim radius();
};

#endif