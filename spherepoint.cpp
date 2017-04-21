#include "spherepoint.h"
SpherePoint::SpherePoint(Loc sphere_loc, Secs secs):Point(secs){
    _sphere_loc = sphere_loc;
}

SpherePoint::~SpherePoint(){}

Loc SpherePoint::sphere_loc(){
    return _sphere_loc;
}