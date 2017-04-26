#include "infcylinder.h"

InfCylinder::InfCylinder(Dim radius, string mat, Kelvin t_init): Geometry(mat,t_init){
	_radius = radius;
	vector<InfCylinderPoint> _temp_dist;	
}

InfCylinder::InfCylinder(Dim radius, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
	_radius = radius;
}

InfCylinder::~InfCylinder(){}

float InfCylinder::radius(){
	return _radius;
}

void InfCylinder::temp_dist(vector<InfCylinderPoint> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<InfCylinderPoint> &InfCylinder::temp_dist() {
	return _temp_dist;
}