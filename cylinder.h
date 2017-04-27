#ifndef CYLINDER_H
#define CYLINDER_H
#include "geometry.h"
#include "cylinderpoint.h"

class Cylinder : public Geometry{
	Dim _radius;
	Dim _length;
	vector<CylinderPoint> _temp_dist;
public:
	Cylinder(Dim radius, Dim length, string mat, Kelvin t_init);
	Cylinder(Dim radius, Dim length, float k, float c, float p, Kelvin t_init);
	~Cylinder();
	Dim radius();
	Dim length();
	void temp_dist(vector<CylinderPoint> &temp_dist); 
	vector<CylinderPoint> &temp_dist(); 
	bool validpoint(CylinderPoint &p);
};

#endif