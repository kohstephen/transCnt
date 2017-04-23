#ifndef ENVMAT_H
#define ENVMAT_H
#include "constant.h"

class EnvMat{
    string _envmat;
    Temp _t_inf;
    float _h;
public:
    EnvMat(string envmat, Temp t_inf);
    ~EnvMat();
    float h();
    float t_inf();
};

#endif