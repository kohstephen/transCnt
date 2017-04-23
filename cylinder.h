#ifndef CYLINDER_H
#define CYLINDER_H
#include "geometry.h"

class Cylinder : public Geometry{
	Dim _radius;
	Dim _length;
public:
	Cylinder(Dim radius, Dim length, string mat, Temp t_init);
	~Cylinder();
	Dim radius();
	Dim length();
};

#endif