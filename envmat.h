#ifndef ENVMAT_H
#define ENVMAT_H
#include "constant.h"

class EnvMat{
    Kelvin _t_inf;
    float _h;
public:
    EnvMat(string envmat, Kelvin t_inf);
    EnvMat(float h, Kelvin t_inf);
    ~EnvMat();
    float h();
    Kelvin t_inf();
};

#endif