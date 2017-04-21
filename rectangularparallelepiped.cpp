#include "rectangularparallelepiped.h"

RectangularParallelepiped::RectangularParallelepiped(float L_1, float L_2, float L_3, string mat, Temp t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
	_L_3 = L_3;
}

RectangularParallelepiped::~RectangularParallelepiped(){}

float RectangularParallelepiped::getL1(){
	return _L_1;
}

float RectangularParallelepiped::getL2(){
	return _L_2;
}

float RectangularParallelepiped::getL3(){
	return _L_3;
}