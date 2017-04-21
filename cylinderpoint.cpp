#include "cylinderpoint.h"

CylinderPoint::CylinderPoint(Loc cyl_loc, Loc rect_loc, Secs secs):Point(secs){
    _cyl_loc = cyl_loc;
    _rect_loc = rect_loc;
}

CylinderPoint::~CylinderPoint(){}

Loc CylinderPoint::cyl_loc(){
    return _cyl_loc;
}

Loc CylinderPoint::rect_loc(){
    return _rect_loc;
}