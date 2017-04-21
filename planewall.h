#ifndef PLANEWALL_H
#define PLANEWALL_H
#include "constant.h"
#include "geometry.h"
#include <string>

class PlaneWall : public Geometry{
	float _length;

public:
	PlaneWall(float length, string mat, Temp t_init);
	~PlaneWall();
	float length();
};

#endif