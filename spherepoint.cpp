#include "spherepoint.h"
SpherePoint::SpherePoint(Loc sphere_loc, Secs secs){
    _sphere_loc = sphere_loc;
    _secs = secs;
    _temp = 0.0f;
}

SpherePoint::~SpherePoint(){}

Loc SpherePoint::sphere_loc(){
    return _sphere_loc;
}

Temp SpherePoint::temp(){
    return _temp;
}

void SpherePoint::temp(Temp temp){
    _temp = temp;
}

Secs SpherePoint::time(){
    return _secs;
}