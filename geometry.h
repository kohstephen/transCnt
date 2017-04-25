#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "constant.h"

class Geometry{
private:
	Kelvin _t_init;
	float _k;
	float _c;
	float _p;
	float _a;
	void calculate(string mat);
	
public:
	Geometry(string mat, Kelvin t_init);
	Geometry(float k, float c, float p, Kelvin t_init);
	~Geometry();
	Kelvin t_init();
	float k();
	float c();
	float p();
	float a();
};

#endif