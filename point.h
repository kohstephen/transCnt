#ifndef POINT_H
#define POINT_H
#include "constant.h"

class Point{
	Temp _temp;
	Secs _secs;
public:
	Point(Secs secs);
	~Point();
	Temp temp();
	void temp(Temp temp);
	Secs time();
};

#endif