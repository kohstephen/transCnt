#ifndef POINT_H
#define POINT_H
#include "constant.h"

class Point{
	Kelvin _temp;
	Secs _secs;
public:
	Point(Secs secs);
	~Point();
	Kelvin temp();
	void temp(Kelvin temp);
	Secs time();
};

#endif