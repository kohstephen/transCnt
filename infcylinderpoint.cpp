#include "infcylinderpoint.h"

InfCylinderPoint::InfCylinderPoint(Loc cyl_loc, Secs secs) : Point(secs){
    _cyl_loc = cyl_loc;
}

InfCylinderPoint::~InfCylinderPoint(){}

Loc InfCylinderPoint::cyl_loc(){
    return _cyl_loc;
}