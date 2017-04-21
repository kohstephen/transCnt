#ifndef RECTBAR_H
#define RECTBAR_H
#include "geometry.h"

class RectBar : public Geometry{
	float _L_1;
	float _L_2;
	float _L_3;
public:
	RectBar(float L_1, float L_2, float L_3, string mat, Temp t_init);
	~RectBar();
	float getL1();
	float getL2();
	float getL3();
};

#endif