#ifndef INFCYLINDER_H
#define INFCYLINDER_H
#include "geometry.h"
#include "infcylinderpoint.h"

class InfCylinder : public Geometry{
	Dim _radius;
	vector<InfCylinderPoint> _temp_dist;
public:
	InfCylinder(Dim radius, string mat, Kelvin t_init);
	InfCylinder(Dim radius, float k, float c, float p, Kelvin t_init);
	~InfCylinder();
	Dim radius();
	void temp_dist(vector<InfCylinderPoint> &temp_dist); 
	vector<InfCylinderPoint> &temp_dist(); 
};

#endif