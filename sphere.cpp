#include "sphere.h"

Sphere::Sphere(Dim radius, string mat, Kelvin t_init) : Geometry(mat,t_init){
		_radius = radius;
		vector<SpherePoint> _temp_dist;	
}

Sphere::Sphere(Dim radius, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
	_radius = radius;
}

Sphere::~Sphere(){}

Dim Sphere::radius(){
	return _radius;
}

void Sphere::temp_dist(vector<SpherePoint> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<SpherePoint> &Sphere::temp_dist() {
	return _temp_dist;
}