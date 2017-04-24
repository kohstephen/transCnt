#ifndef PLANEWALL_H
#define PLANEWALL_H
#include "constant.h"
#include "geometry.h"
#include <string>
#include "planewallpoint.h"

class PlaneWall : public Geometry{
	Dim _length;
	vector<PlaneWallPoint> _temp_dist;

public:
	PlaneWall(Dim length, string mat, Kelvin t_init);
	~PlaneWall();
	Dim length();
	void temp_dist(vector<PlaneWallPoint> &temp_dist); 
};

#endif