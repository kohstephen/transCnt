#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "constant.h"

class Geometry{
private:
	string _mat;
	Temp _t_init;
	float _k;
	float _c;
	float _p;
	float _a;
	void calculate();
	
public:
	Geometry(string mat, Temp t_init);
	~Geometry();
	string mat();
	Temp t_init();
	float k();
	float c();
	float p();
	float a();
};

#endif