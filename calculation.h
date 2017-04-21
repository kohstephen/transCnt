#ifndef CALCULATION_H
#define CALCULATION_H
#include <vector>
#include <string>
#include "planewall.h"
#include "infinitecylinder.h"
#include "rectangularparallelepiped.h"
#include "cylinder.h"
#include "infiniterectangularbar.h"
#include "sphere.h"
#include "constant.h"
#include "planewallpoint.h"
#include "spherepoint.h"

using namespace std;

extern const float PI;

/**
void temp_at_time(vector<float>* ret, Sphere s, string mat, string envmat, 
    vector<float>& points, float time, float t_init, float t_inf);
**/
void temp_at_point(PlaneWall &w, PlaneWallPoint &p, string envmat, Temp t_inf);

void temp_at_point(Sphere &s, SpherePoint &p, string envmat, Temp t_inf);

/*


float temp_at_time_at_point(InfiniteCylinder cylinder, string mat, string envmat, 
    float x, float time, float t_init, float t_inf);

float temp_at_time_at_point(RectangularParallelepiped rp, string mat, string envmat, 
    float x1, float x2, float x3, float time, float t_init, float t_inf);

float temp_at_time_at_point(Cylinder cylinder, string mat, string envmat, 
    float r, float x, float time, float t_init, float t_inf);

float temp_at_time_at_point(InfiniteRectangularBar irb, string mat, string envmat, 
    float x1, float x2, float time, float t_init, float t_inf);

    
    */
    
#endif