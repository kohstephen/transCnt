#include "geometry.h"
Geometry::Geometry(string mat, Temp t_init){
    _mat = mat;
    _t_init = t_init;
    calculate();
}

Geometry::~Geometry(){}

string Geometry::mat(){
    return _mat;
}

Temp Geometry::t_init(){
    return _t_init;
}

void Geometry::calculate(){
    _k = get_k(_mat, _t_init);
    _c = get_c(_mat, _t_init);
    _p = get_p(_mat);
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