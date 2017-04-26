#ifndef POINT_H
#define POINT_H
#include "constant.h"

/**
 * Point class that represent the basics of a point.
 * Other specific point classes are subclasses of Point.
 */
class Point{
	//! temparature at the point; user does not set this
	Kelvin _temp;
	//! Every point should have a special dimension: time
	Secs _secs;
public:
	Point(Secs secs);
	~Point();
	Kelvin temp();
	void temp(Kelvin temp);
	Secs time();
};

#endif