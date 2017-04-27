#include "cylinder.h"

Cylinder::Cylinder(Dim radius, Dim length, string mat, Kelvin t_init):Geometry(mat,t_init){
	_radius = radius;
	_length = length;
	vector<CylinderPoint> _temp_dist;	
}

Cylinder::Cylinder(Dim radius, Dim length, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
	_radius = radius;
	_length = length;
}


Cylinder::~Cylinder(){}

Dim Cylinder::radius(){
	return _radius;
}

Dim Cylinder::length(){
	return _length;
}

void Cylinder::temp_dist(vector<CylinderPoint> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<CylinderPoint> &Cylinder::temp_dist() {
	return _temp_dist;
}

bool Cylinder::validpoint(CylinderPoint &p){
	Loc cyll = p.cyl_loc();
	Loc rectl = p.rect_loc();
	if(cyll<0 || cyll>_radius || rectl<0 || rectl>_length) return false;
	return true;
}
