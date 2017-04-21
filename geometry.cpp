#include "geometry.h"
Geometry::Geometry(string mat, Temp t_init){
    _mat = mat;
    _t_init = t_init;
}

Geometry::~Geometry(){}

string Geometry::mat(){
    return _mat;
}

Temp Geometry::t_init(){
    return _t_init;
}