#ifndef SPHERE_H
#define SPHERE_H
#include "geometry.h"
#include "spherepoint.h"

class Sphere : public Geometry{
	Dim _radius;
	vector<SpherePoint> _temp_dist;	

public:
	Sphere(Dim radius, string mat, Kelvin t_init);
	Sphere(Dim radius, float k, float c, float p, Kelvin t_init);
	~Sphere();
	Dim radius();
	void temp_dist(vector<SpherePoint> &temp_dist); 
	vector<SpherePoint> &temp_dist(); 
	//! check if the SpherePoint p is a valid point inside this Sphere
	bool validpoint(SpherePoint &p);
};


#endif