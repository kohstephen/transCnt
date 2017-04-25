#include "envmat.h"

EnvMat::EnvMat(string envmat, Kelvin t_inf){
    _t_inf = t_inf;
    _h = get_h(envmat);
}

EnvMat::EnvMat(float h, Kelvin t_inf){
    _t_inf = t_inf;
    _h = h;
}

EnvMat::~EnvMat(){}

float EnvMat::h(){
    return _h;
}

Kelvin EnvMat::t_inf(){
    return _t_inf;
}