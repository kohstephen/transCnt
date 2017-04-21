#include "planewall.h"
#include <vector>
#include <utility>

// mesh_density is the number of regions 
PlaneWall::PlaneWall(float length, string mat, Temp t_init) : Geometry(mat,t_init){ // float mesh_density = 0){
	_length = length;

	// extra 1 to cover both ends of divided region
	// Vector<Point> temp_dist (mesh_density + 1) = ;	
}

PlaneWall::~PlaneWall(){}

float PlaneWall::length(){
	return _length;
}