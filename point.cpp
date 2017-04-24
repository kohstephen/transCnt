#include "point.h"

Point::Point(Secs secs){
    _secs = secs;
    _temp = 0.0f;
}

Point::~Point(){}

Kelvin Point::temp(){
    return _temp;
}

void Point::temp(Kelvin temp){
    _temp = temp;
}

Secs Point::time(){
    return _secs;
}