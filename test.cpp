#include <math.h>
#include <iostream>
#include "calculation.h"
#include <chrono>

int main(){
    string mat = "st";
    Kelvin t_init = 500;

    EnvMat envmat = EnvMat("water", 300);
    

    //(zeta = 1.142)
    
    // Sphere s = Sphere(0.025);
    // output = temp_at_time_at_point(s, mat, envmat, r, time, t_init, t_inf);
    
    // cout << output << endl;
    // //assert( abs(output-result) < 0.01 );

    Sphere s2 = Sphere(0.001, mat, t_init);
    SpherePoint sp2 = SpherePoint(0.0005,5);
    temp_at_point(s2, sp2, envmat);
    //cout << sp2.temp() << endl;

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

    InfCylinder ic3 = InfCylinder(0.1,mat,t_init);
    InfCylinderPoint icp3 = InfCylinderPoint(0.05,100);
    temp_at_point(ic3, icp3, envmat);
    //cout << icp3.temp() << endl;


    PlaneWall w = PlaneWall(0.001, mat, t_init);
    PlaneWallPoint wp = PlaneWallPoint(0.0005f, 10.0f);
    temp_at_point(w, wp, envmat);
    temp_on_mesh(w, 10, 10, envmat); 
    cout << endl << wp.temp() << endl << endl;

    PlaneWall w2 = PlaneWall(0.01, mat, t_init);
    PlaneWallPoint wp2 = PlaneWallPoint(0.005f, 10.0f);
    temp_at_point(w2, wp2, envmat);
    temp_on_mesh(w2, 10, 10, envmat); 
    cout << endl << wp2.temp() << endl << endl;

    PlaneWall w3 = PlaneWall(0.1, mat, t_init);
    PlaneWallPoint wp3 = PlaneWallPoint(0.05f, 100.0f);
    temp_at_point(w3, wp3, envmat);
    temp_on_mesh(w3, 100, 10, envmat); 
    cout << endl << wp3.temp() << endl << endl;

    PlaneWall w4 = PlaneWall(0.1, mat, t_init);
    PlaneWallPoint wp4 = PlaneWallPoint(0.05f, 10.0f);
    temp_at_point(w4, wp4, envmat);
    temp_on_mesh(w4, 10, 10, envmat);
    cout << endl << wp4.temp() << endl << endl;
    
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