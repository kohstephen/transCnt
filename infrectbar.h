#ifndef INFRECTBAR_H
#define INFRECTBAR_H
#include "geometry.h"

class InfRectBar : public Geometry{
	float _L_1;
	float _L_2;
public:
	InfRectBar(float L_1, float L_2, string mat, Temp t_init);
	~InfRectBar();
	float getL1();
	float getL2();
};

#endif