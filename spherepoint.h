#ifndef SPHEREPOINT_H
#define SPHEREPOINT_H
#include "point.h"

class SpherePoint : public Point{
private:
	Loc _sphere_loc;

public:
	SpherePoint(Loc sphere_loc, Secs secs);
	~SpherePoint();
	Loc sphere_loc();
};

#endif