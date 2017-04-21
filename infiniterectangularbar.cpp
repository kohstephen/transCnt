#include "infiniterectangularbar.h"

InfiniteRectangularBar::InfiniteRectangularBar(float L_1, float L_2, string mat, Temp t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
}

InfiniteRectangularBar::~InfiniteRectangularBar(){}

float InfiniteRectangularBar::getL1(){
	return _L_1;
}

float InfiniteRectangularBar::getL2(){
	return _L_2;
}
