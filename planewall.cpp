#include "planewall.h"
#include <vector>
#include <utility>
#include "point.h"
#include <iostream>

// mesh_density is the number of regions 
PlaneWall::PlaneWall(Dim length, string mat, Kelvin t_init) : Geometry(mat,t_init){ 
	_length = length;
	vector<PlaneWallPoint> _temp_dist;	
}

PlaneWall::PlaneWall(Dim length, float k, float c, float p, Kelvin t_init):Geometry(k,c,p,t_init){
	_length = length;
}

PlaneWall::~PlaneWall(){}

Dim PlaneWall::length(){
	return _length;
}

void PlaneWall::temp_dist(vector<PlaneWallPoint> &temp_dist) {
   _temp_dist = temp_dist; 
}
