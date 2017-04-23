#ifndef INFRECTBAR_H
#define INFRECTBAR_H
#include "geometry.h"

class InfRectBar : public Geometry{
	Dim _L_1;
	Dim _L_2;
public:
	InfRectBar(Dim L_1, Dim L_2, string mat, Temp t_init);
	~InfRectBar();
	Dim l1();
	Dim l2();
};

#endif