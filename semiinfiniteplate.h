#ifndef SEMIINFINITEPLATE_H
#define SEMIINFINITEPLATE_H
#include "geometry.h"

class SemiInfinitePlate : public Geometry{
	float _L;
public:
	SemiInfinitePlate(float L, string mat, Temp t_init);
	~SemiInfinitePlate();
	float getL();
};

#endif