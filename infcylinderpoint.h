#ifndef INFCYLINDERPOINT_H
#define INFCYLINDERPOINT_H
#include "point.h"

class InfCylinderPoint : public Point{
	Loc _cyl_loc;

public:
	InfCylinderPoint(Loc cyl_loc, Secs secs);
	~InfCylinderPoint();
	Loc cyl_loc();
};

#endif