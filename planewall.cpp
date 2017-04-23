#include "planewall.h"
#include <vector>
#include <utility>
#include "point.h"
#include <iostream>

// mesh_density is the number of regions 
PlaneWall::PlaneWall(Dim length, string mat, Temp t_init) : Geometry(mat,t_init){ 
	_length = length;
	vector<PlaneWallPoint> _temp_dist;	
}

PlaneWall::~PlaneWall(){}

Dim PlaneWall::length(){
	return _length;
}

/*
friend void PlaneWall::temp_dist(int num_points) {
    float incr = _length / (num_points - 1.0f); 
    int counter = 0;
    vector<PlaneWallPoint> temp_dist (num_points);
    for (auto it = temp_dist.begin(); it != temp_dist.end(); it++) {
    	(*it).rect_loc(counter*incr);	
    	counter++;
    }
    // to check rect_loc are correctly spaced
    for (auto i: temp_dist)
        cout << i.rect_loc() << ' ';
    _temp_dist = temp_dist; 
}
*/