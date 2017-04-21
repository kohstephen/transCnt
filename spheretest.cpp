#include <math.h>
#include <iostream>
#include "calculation.h"
#include <chrono>

int main(){
    string mat = "st";
    string envmat = "water";
    Temp t_init = 500;
    Temp t_inf = 300;
    

    //(zeta = 1.142)
    
    // Sphere s = Sphere(0.025);
    // output = temp_at_time_at_point(s, mat, envmat, r, time, t_init, t_inf);
    
    // cout << output << endl;
    // //assert( abs(output-result) < 0.01 );

    // Sphere s2 = Sphere(0.001);
    // output = temp_at_time_at_point(s2, mat, envmat, 0.0005, 5, t_init, t_inf);
    // cout << output << endl;

    // Sphere s3 = Sphere(0.01);
    // output = temp_at_time_at_point(s3, mat, envmat, 0.005, 10, t_init, t_inf);
    // cout << output << endl;

    // Sphere s4 = Sphere(0.1);
    // output = temp_at_time_at_point(s4, mat, envmat, 0.05, 100, t_init, t_inf);
    // cout << output << endl;

    // InfiniteCylinder ic1 = InfiniteCylinder(0.001);
    // output = temp_at_time_at_point(ic1, mat, envmat, 0.0005, 10, t_init, t_inf);
    // cout << output << endl;

    // InfiniteCylinder ic2 = InfiniteCylinder(0.01);
    // output = temp_at_time_at_point(ic2, mat, envmat, 0.005, 10, t_init, t_inf);
    // cout << output << endl;

    // InfiniteCylinder ic3 = InfiniteCylinder(0.1);
    // output = temp_at_time_at_point(ic3, mat, envmat, 0.05, 100, t_init, t_inf);
    // cout << output << endl;


    PlaneWall w = PlaneWall(0.001, mat, t_init);
    PlaneWallPoint_impl p = PlaneWallPoint_impl(0.0005f, 10.0f);
    temp_at_point(w, p, envmat, t_inf);
    cout << p.temp() << endl;
    //cout << p._temp << endl;


    // PlaneWall pl2 = PlaneWall(0.01);
    // output = temp_at_time_at_point(pl2, mat, envmat, 0.005, 10, t_init, t_inf);
    // cout << output << endl;

    // PlaneWall pl3 = PlaneWall(0.1);
    // output = temp_at_time_at_point(pl3, mat, envmat, 0.05, 100, t_init, t_inf);
    // cout << output << endl;

    // InfiniteCylinder ic = InfiniteCylinder(0.001);
    // output = temp_at_time_at_point(ic, mat, envmat, 0.0005, 10, t_init, t_inf);
    // cout << output << endl;


}
