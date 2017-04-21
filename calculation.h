#ifndef CALCULATION_H
#define CALCULATION_H
#include <vector>
#include <string>
#include "planewall.h"
#include "planewallpoint.h"
#include "infcylinder.h"
#include "infcylinderpoint.h"
#include "rectbar.h"
#include "rectbarpoint.h"
#include "cylinder.h"
#include "cylinderpoint.h"
#include "infrectbar.h"
#include "infrectbarpoint.h"
#include "sphere.h"
#include "spherepoint.h"
#include "constant.h"

using namespace std;

extern const float PI;

/**
void temp_at_time(vector<float>* ret, Sphere s, string mat, string envmat, 
    vector<float>& points, float time, float t_init, float t_inf);
**/

void temp_at_point(PlaneWall &w, PlaneWallPoint &p, string envmat, Temp t_inf);

void temp_at_point(Sphere &s, SpherePoint &p, string envmat, Temp t_inf);

void temp_at_point(InfCylinder &icyl, InfCylinderPoint &p, string envmat, Temp t_inf);

void temp_at_point(RectBar &rb, RectBarPoint &p, string envmat, Temp t_inf);

void temp_at_point(Cylinder &cyl, CylinderPoint &p, string envmat, Temp t_inf);

void temp_at_point(InfRectBar &irb, InfRectBarPoint &p, string envmat, Temp t_inf);

    
#endif