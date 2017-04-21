#include "semiinfiniteplate.h"

SemiInfinitePlate::SemiInfinitePlate(float L, string mat, Temp t_init):Geometry(mat,t_init){
	_L = L;
}
SemiInfinitePlate::~SemiInfinitePlate(){}

float SemiInfinitePlate::getL(){
	return _L;
}