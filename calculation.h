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
#include "envmat.h"

using namespace std;

extern const float PI;

/**
 * Calculate temperature at a certain point.
 * The temperature will be set in the point provided by the user.
 */
void temp_at_point(PlaneWall &w, PlaneWallPoint &p, EnvMat &envmat);

void temp_at_point(Sphere &s, SpherePoint &p, EnvMat &envmat);

void temp_at_point(InfCylinder &icyl, InfCylinderPoint &p, EnvMat &envmat);

void temp_at_point(RectBar &rb, RectBarPoint &p, EnvMat &envmat);

void temp_at_point(Cylinder &cyl, CylinderPoint &p, EnvMat &envmat);

void temp_at_point(InfRectBar &irb, InfRectBarPoint &p, EnvMat &envmat);


void temp_on_mesh(PlaneWall &w, Secs secs, int mesh_density, EnvMat &envmat); 

void temp_on_mesh(InfCylinder &icyl, Secs secs, int mesh_density, EnvMat &envmat); 

void temp_on_mesh(Sphere &s, Secs secs, int mesh_density, EnvMat &envmat); 

void temp_on_mesh(InfRectBar &irb, Secs secs, int mesh_density, EnvMat &envmat);

void temp_on_mesh(Cylinder &cyl, Secs secs, int mesh_density, EnvMat &envmat);

void temp_on_mesh(RectBar &rb, Secs secs, int mesh_density, EnvMat &envmat);

/**
 * Utility functions to convert Kelvin to Fahrenheit and Celcius,
 * and from Fahrenheit and Celcius to Kelvin.
 */
float kelvin_to_fahrenheit(Kelvin k);

float kelvin_to_celcius(Kelvin k);

Kelvin  fahrenheit_to_kelvin(float f);

Kelvin celcius_to_kelvin(float c);

// avg
float avg_temp_at_time(Sphere &s, Secs time, EnvMat &envmat);

float avg_temp_at_time(PlaneWall &s, Secs time, EnvMat &envmat);

float avg_temp_at_time(InfCylinder &s, Secs time, EnvMat &envmat);

float avg_temp_at_time(Cylinder &s, Secs time, EnvMat &envmat);

float avg_temp_at_time(RectBar &s, Secs time, EnvMat &envmat);

float avg_temp_at_time(InfRectBar &s, Secs time, EnvMat &envmat);

// Plot
void plot(Sphere &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(PlaneWall &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(InfCylinder &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(Cylinder &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(RectBar &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(InfRectBar &s, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(Sphere &s, vector<SpherePoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(PlaneWall &s, vector<PlaneWallPoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(InfCylinder &s, vector<InfCylinderPoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(Cylinder &s, vector<CylinderPoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(RectBar &s, vector<RectBarPoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);

void plot(InfRectBar &s, vector<InfRectBarPoint> &p, Secs start, Secs end, Secs intrv, EnvMat &envmat);
#endif
