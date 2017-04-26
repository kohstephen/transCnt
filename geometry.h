#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "constant.h"
// #include "point.h"

/**
 * Geometry class that represent the basics of a geometry.
 * Other specific geometry classes are subclasses of Geometry.
 */
class Geometry{
private:
	//! initial temparature of the geometry
	Kelvin _t_init;
	//! k,c,p,a value stored
	float _k;
	float _c;
	float _p;
	float _a;
	//! private function to calculate k,c,p,a value from material string
	void calculate(string mat);
	// vector<Point> _temp_dist;	
	
public:
	/**
	 * User should provide the initial temparature,
	 * plus either a material name,
	 * or the k,c,p values if they want a customized material
	 * when instantiating a geometry.
	 */
	Geometry(string mat, Kelvin t_init);
	Geometry(float k, float c, float p, Kelvin t_init);
	~Geometry();
	Kelvin t_init();
	float k();
	float c();
	float p();
	float a();
	// void temp_dist(vector<Point> &temp_dist); 
	// vector<Point> &temp_dist(); 
};

#endif