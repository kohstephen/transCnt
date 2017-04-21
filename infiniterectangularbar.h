#ifndef INFINITERECTANGULARBAR_H
#define INFINITERECTANGULARBAR_H
#include "geometry.h"

class InfiniteRectangularBar : public Geometry{
	float _L_1;
	float _L_2;
public:
	InfiniteRectangularBar(float L_1, float L_2, string mat, Temp t_init);
	~InfiniteRectangularBar();
	float getL1();
	float getL2();
};

#endif