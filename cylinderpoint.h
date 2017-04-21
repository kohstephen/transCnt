#ifndef CYLINDERPOINT_H
#define CYLINDERPOINT_H
#include "point.h"

class CylinderPoint : public Point{
	Loc _cyl_loc;
	Loc _rect_loc;

public:
	CylinderPoint(Loc cyl_loc, Loc rect_loc, Secs secs);
	~CylinderPoint();
	Loc cyl_loc();
	Loc rect_loc();
};

#endif