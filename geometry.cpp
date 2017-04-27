#include "geometry.h"
Geometry::Geometry(string mat, Kelvin t_init){
    _t_init = t_init;
    calculate(mat);
	// vector<Point> _temp_dist;	
    _method = -1;
}

Geometry::Geometry(float k, float c, float p, Kelvin t_init){
    _t_init = t_init;
    _k = k;
    _c = c;
    _p = p;
    _a =  _k/(_p*_c);
    // vector<Point> _temp_dist;
    _method = -1;
}

Geometry::~Geometry(){}

Kelvin Geometry::t_init(){
    return _t_init;
}

void Geometry::method(int method){
    _method = method;
}

int Geometry::method(){
    return _method;
}

void Geometry::calculate(string mat){
    _k = get_k(mat, _t_init);
    _c = get_c(mat, _t_init);
    _p = get_p(mat);
    _a =  _k/(_p*_c);
}

float Geometry::k(){
    return _k;
}

float Geometry::c(){
    return _c;
}

float Geometry::p(){
    return _p;
}

float Geometry::a(){
    return _a;
}
 /*
void Geometry::temp_dist(vector<Point> &temp_dist) {
   _temp_dist = temp_dist; 
}

vector<Point> &Geometry::temp_dist() {
	return _temp_dist;
}
*/
