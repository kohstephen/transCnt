#include "envmat.h"

EnvMat::EnvMat(string envmat, Temp t_inf){
    _envmat = envmat;
    _t_inf = t_inf;
    _h = get_h(envmat);
}

EnvMat::~EnvMat(){}

float EnvMat::h(){
    return _h;
}

float EnvMat::t_inf(){
    return _t_inf;
}