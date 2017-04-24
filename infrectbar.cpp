#include "infrectbar.h"

InfRectBar::InfRectBar(Dim L_1, Dim L_2, string mat, Kelvin t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
}

InfRectBar::~InfRectBar(){}

Dim InfRectBar::l1(){
	return _L_1;
}

Dim InfRectBar::l2(){
	return _L_2;
}
