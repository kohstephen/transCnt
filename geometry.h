#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "constant.h"

class Geometry{
	string _mat;
	Temp _t_init;
	
public:
	Geometry(string mat, Temp t_init);
	~Geometry();
	string mat();
	Temp t_init();
};

#endif