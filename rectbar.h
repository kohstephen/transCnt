#ifndef RECTBAR_H
#define RECTBAR_H
#include "geometry.h"

class RectBar : public Geometry{
	Dim _L_1;
	Dim _L_2;
	Dim _L_3;
public:
	RectBar(Dim L_1, Dim L_2, Dim L_3, string mat, Kelvin t_init);
	~RectBar();
	Dim l1();
	Dim l2();
	Dim l3();
};

#endif