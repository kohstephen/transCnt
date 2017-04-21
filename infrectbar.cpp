#include "infrectbar.h"

InfRectBar::InfRectBar(float L_1, float L_2, string mat, Temp t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
}

InfRectBar::~InfRectBar(){}

float InfRectBar::getL1(){
	return _L_1;
}

float InfRectBar::getL2(){
	return _L_2;
}
