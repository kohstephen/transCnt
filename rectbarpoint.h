#ifndef RECTBARPOINT_H
#define RECTBARPOINT_H
#include "point.h"

class RectBarPoint : public Point{
private:
	Loc _rect_loc1, _rect_loc2, _rect_loc3;

public:
	RectBarPoint(Loc rect_loc1, rect_loc2, rect_loc3, Secs secs);
	~RectBarPoint();
	Loc rect_loc1();
	Loc rect_loc2();
	Loc rect_loc3();
};

#endif