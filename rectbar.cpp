#include "rectbar.h"

RectBar::RectBar(Dim L_1, Dim L_2, Dim L_3, string mat, Kelvin t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
	_L_3 = L_3;
}

RectBar::~RectBar(){}

Dim RectBar::l1(){
	return _L_1;
}

Dim RectBar::l2(){
	return _L_2;
}

Dim RectBar::l3(){
	return _L_3;
}