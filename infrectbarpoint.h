#ifndef INFRECTBARPOINT_H
#define INFRECTBARPOINT_H
#include "point.h"

/**
 * InfRectBarPoint is a subclass of Point.
 * It represents a point inside infinite rectangular bar.
 */
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