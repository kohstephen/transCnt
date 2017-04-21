#include "rectbarpoint.h"

RectBarPoint::RectBarPoint(Loc rect_loc1, Loc rect_loc2, Loc rect_loc3, Secs secs):Point(secs){
    _rect_loc1 = rect_loc1;
    _rect_loc2 = rect_loc2;
    _rect_loc3 = rect_loc3;
}

RectBarPoint::~RectBarPoint(){}

Loc RectBarPoint::rect_loc1(){
    return _rect_loc1;
}
	
Loc RectBarPoint::rect_loc2(){
    return _rect_loc2;
}

Loc RectBarPoint::rect_loc3(){
    return _rect_loc3;
}