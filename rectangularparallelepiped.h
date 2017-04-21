#ifndef RECTANGULARPARALLELEPIPED_H
#define RECTANGULARPARALLELEPIPED_H
#include "geometry.h"

class RectangularParallelepiped : public Geometry{
	float _L_1;
	float _L_2;
	float _L_3;
public:
	RectangularParallelepiped(float L_1, float L_2, float L_3, string mat, Temp t_init);
	~RectangularParallelepiped();
	float getL1();
	float getL2();
	float getL3();
};

#endif