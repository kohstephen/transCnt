#include "infrectbarpoint.h"

InfRectBarPoint::InfRectBarPoint(Loc rect_loc1, Loc rect_loc2, Secs secs):Point(secs){
    _rect_loc1 = rect_loc1;
    _rect_loc2 = rect_loc2;
}

InfRectBarPoint::~InfRectBarPoint(){}

Loc InfRectBarPoint::rect_loc1(){
    return _rect_loc1;
}

Loc InfRectBarPoint::rect_loc2(){
    return _rect_loc2;
}