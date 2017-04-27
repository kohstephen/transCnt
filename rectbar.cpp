#include "rectbar.h"

RectBar::RectBar(Dim L_1, Dim L_2, Dim L_3, string mat, Kelvin t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
	_L_3 = L_3;
	vector<RectBarPoint> _temp_dist;
}

RectBar::RectBar(Dim L_1, Dim L_2, Dim L_3, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
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

bool RectBar::validpoint(RectBarPoint &p){
	Loc l1 = p.rect_loc1();
	Loc l2 = p.rect_loc2();
	Loc l3 = p.rect_loc3();
	if(l1>_L_1 || l1<0 || l2>_L_2 || l2<0 || l3>_L_3|| l3<0) return false;
	return true;
}

void RectBar::temp_dist(vector<RectBarPoint> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<RectBarPoint> &RectBar::temp_dist() {
	return _temp_dist;
}