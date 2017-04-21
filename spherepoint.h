#ifndef SPHEREPOINT_H
#define SPHEREPOINT_H
#include "constant.h"

class SpherePoint{
private:
	Loc _sphere_loc;
	Secs _secs;
	Temp _temp;
public:
	SpherePoint(Loc sphere_loc, Secs secs);
	~SpherePoint();
	Loc sphere_loc();
	Temp temp();
	void temp(Temp temp);
	Secs time();
};

#endif