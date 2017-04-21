#ifndef INFRECTBARPOINT_H
#define INFRECTBARPOINT_H
#include "point.h"

class InfRectBarPoint : public Point{
	Loc _rect_loc1;
	Loc _rect_loc2;

public:
	InfRectBarPoint(Loc rect_loc1, Loc rect_loc2, Secs secs);
	~InfRectBarPoint();
	Loc rect_loc1();
	Loc rect_loc2();
};

#endif