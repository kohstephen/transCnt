#include "point.h"

Point::Point(Secs secs){
    _secs = secs;
    _temp = 0.0f;
}

Point::~Point(){}

Temp Point::temp(){
    return _temp;
}

void Point::temp(Temp temp){
    _temp = temp;
}

Secs Point::time(){
    return _secs;
}