#include "rectangularparallelepiped.h"

RectBar::RectBar(float L_1, float L_2, float L_3, string mat, Temp t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
	_L_3 = L_3;
}

RectBar::~RectBar(){}

float RectBar::getL1(){
	return _L_1;
}

float RectBar::getL2(){
	return _L_2;
}

float RectBar::getL3(){
	return _L_3;
}