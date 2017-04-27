#ifndef RECTBAR_H
#define RECTBAR_H
#include "geometry.h"
#include "rectbarpoint.h"

class RectBar : public Geometry{
	Dim _L_1;
	Dim _L_2;
	Dim _L_3;
	vector<RectBarPoint> _temp_dist;
public:
	RectBar(Dim L_1, Dim L_2, Dim L_3, string mat, Kelvin t_init);
	RectBar(Dim L_1, Dim L_2, Dim L_3, float k, float c, float p, Kelvin t_init);
	~RectBar();
	Dim l1();
	Dim l2();
	Dim l3();
	//! check if the RectBarPoint is a valid point for this geometry
	bool validpoint(RectBarPoint &p);
	void temp_dist(vector<RectBarPoint> &temp_dist); 
	vector<RectBarPoint> &temp_dist(); 
};

#endif