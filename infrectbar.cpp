#include "infrectbar.h"

InfRectBar::InfRectBar(Dim L_1, Dim L_2, string mat, Kelvin t_init):Geometry(mat,t_init){
	_L_1 = L_1;
	_L_2 = L_2;
	vector<InfRectBarPoint> _temp_dist;	
}

InfRectBar::InfRectBar(Dim L_1, Dim L_2, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
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

void InfRectBar::temp_dist(vector<InfRectBarPoint> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<InfRectBarPoint> &InfRectBar::temp_dist() {
	return _temp_dist;
}

bool InfRectBar::validpoint(InfRectBarPoint &p){
	Loc l1 = p.rect_loc1();
	Loc l2 = p.rect_loc2();
	if(l1>_L_1 || l1<0 || l2>_L_2 || l2<0) return false;
	return true;
}