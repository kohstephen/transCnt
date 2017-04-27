#include <math.h>
#include <iostream>
#include "calculation.h"
#include <chrono>
#include "assert.h"

/*
   T_0 = 500 K steel in T_inf = 300 K water

   LUMPED CAP
   plane wall, L = 0.001, T(10) = 317.529 K
   inf cylinder, r_0 = 0.001, T(10) = 301.536 K
   sphere, r_0 = 0.001, T(10) = T(5) = 305.189 K

   ONE-TERM APPROX
   plane wall, L = 0.01, x = 0.005, T(10) = 460.188 K
   inf cylinder, r_0 = 0.01, r = 0.005, T(10) = 428.534 K
   sphere, r_0 = 0.01, r = 0.005, T(10) = 402.222 K

   MULTI-TERM APPROX 
   plane wall, L = 0.1, x = 0.05, T(100) = 477.353 K (3-term)
   inf cylinder, r_0 = 0.1, r = 0.05, T(100) = 463.285 K (3-term)
   sphere, r_0 = 0.1, r = 0.05, T(100) = 446.296 K (3-term)

   SEMI-INFINITE
   plane wall, L = 0.1, x = 0.05, T(10) = 309.848 K
   inf cylinder, r_0 = 0.1, x = 0.05, T(10) = 309.848 K 
   sphere, r_0 = 0.1, x = 0.05, T(10) = 309.848 K
 */

int main(){
    float ROUNDING = 0.1;
    string mat = "st";
    Kelvin t_init = 500;
    EnvMat envmat = EnvMat("water", 300);

    /*
    InfRectBar irb = InfRectBar(1, 2, mat, t_init);
    temp_on_mesh(irb, 5, 10, envmat); 

    Cylinder cyl = Cylinder(1, 2, mat, t_init);
    temp_on_mesh(cyl, 5, 10, envmat); 
    */ 
    
    RectBar rb = RectBar(.09, .1, .11, mat, t_init);
    vector<int> meshes {9, 20, 45, 99, 214};
    vector<string> ann {"1,000", "10,000", "100,000", "1,000,000", "10M"};

    typedef std::chrono::high_resolution_clock Clock; 

    for (int i = 0; i < meshes.size(); i++) { 
        auto t1 = Clock::now();
        temp_on_mesh(rb, 100, meshes[i], envmat); 
        auto t2 = Clock::now();
        cout << ann[i] + " points: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " ms!" << endl;  
    }

    /**
     * Lumped Capacitance Unit Tests
     */

    Sphere s1 = Sphere(0.001, mat, t_init);
    SpherePoint sp1 = SpherePoint(0.0005,5);
    temp_at_point(s1, sp1, envmat);
    assert(abs(sp1.temp() - 305.189)<ROUNDING);
    temp_on_mesh(s1, 5, 10, envmat); 
    assert(abs(sp1.temp() - s1.temp_dist()[5].temp()) < ROUNDING);    

    InfCylinder ic1 = InfCylinder(0.001, mat, t_init);
    InfCylinderPoint icp1 = InfCylinderPoint(0.0005, 10);
    temp_at_point(ic1, icp1, envmat);
    assert(abs(icp1.temp() - 301.536)<ROUNDING);
    temp_on_mesh(ic1, 10, 10, envmat); 
    assert(abs(icp1.temp() - ic1.temp_dist()[5].temp()) < ROUNDING);    

    PlaneWall w1 = PlaneWall(0.001, mat, t_init);
    PlaneWallPoint wp1 = PlaneWallPoint(0.0005f, 10.0f);
    temp_at_point(w1, wp1, envmat);
    assert(abs(wp1.temp()-317.529)<ROUNDING);
    temp_on_mesh(w1, 10, 10, envmat); 
    assert(abs(wp1.temp() - w1.temp_dist()[5].temp()) < ROUNDING);    


    /**
     * One-Term Approximation Unit Tests
     */

    Sphere s2 = Sphere(0.01, mat, t_init);
    SpherePoint sp2 = SpherePoint(0.005,10);
    temp_at_point(s2, sp2, envmat);
    assert(abs(sp2.temp()-402.222)<ROUNDING);
    temp_on_mesh(s2, 10, 10, envmat); 
    assert(abs(sp2.temp() - s2.temp_dist()[5].temp()) < ROUNDING);    

    // Plot
    //plot(s2, 0, 20, .1, envmat);

    InfCylinder ic2 = InfCylinder(0.01, mat, t_init);
    InfCylinderPoint icp2 = InfCylinderPoint(0.005,10);
    temp_at_point(ic2, icp2, envmat);
    assert(abs(icp2.temp()-428.534)<ROUNDING);
    temp_on_mesh(ic2, 10, 10, envmat); 
    assert(abs(icp2.temp() - ic2.temp_dist()[5].temp()) < ROUNDING);    

    PlaneWall w2 = PlaneWall(0.01, mat, t_init);
    PlaneWallPoint wp2 = PlaneWallPoint(0.005f, 10.0f);
    temp_at_point(w2, wp2, envmat);
    assert(abs(wp2.temp()-460.188)<ROUNDING);
    temp_on_mesh(w2, 10, 10, envmat);
    assert(abs(wp2.temp() - w2.temp_dist()[5].temp()) < ROUNDING);    


    /**
     * Multiple-Term Approximation Unit Tests
     */

    Sphere s3 = Sphere(0.1, mat, t_init);
    SpherePoint sp3 = SpherePoint(0.05,100);
    temp_at_point(s3, sp3, envmat);
    assert(abs(sp3.temp()-446.296)<ROUNDING);
    temp_on_mesh(s3, 100, 10, envmat); 
    assert(abs(sp3.temp() - s3.temp_dist()[5].temp()) < ROUNDING);   

    InfCylinder ic3 = InfCylinder(0.1,mat,t_init);
    InfCylinderPoint icp3 = InfCylinderPoint(0.05,100);
    temp_at_point(ic3, icp3, envmat);
    assert(abs(icp3.temp()-463.285)<ROUNDING);
    temp_on_mesh(ic3, 100, 10, envmat); 
    assert(abs(icp3.temp() - ic3.temp_dist()[5].temp()) < ROUNDING);    
    // Plot
    //plot(ic3, 0, 100, 1, envmat);

    PlaneWall w3 = PlaneWall(0.1, mat, t_init);
    PlaneWallPoint wp3 = PlaneWallPoint(0.05f, 100.0f);
    temp_at_point(w3, wp3, envmat);
    assert(abs(wp3.temp()-477.353)<ROUNDING);
    temp_on_mesh(w3, 100, 10, envmat); 
    assert(abs(wp3.temp() - w3.temp_dist()[5].temp()) < ROUNDING);  

    /**
     * Semi-Infinite Approximation Unit Tests
     */

    Sphere s4 = Sphere(0.1, mat, t_init);
    SpherePoint sp4 = SpherePoint(0.05,10);
    temp_at_point(s4, sp4, envmat);
    assert(abs(sp4.temp()-499.973)<ROUNDING);
    temp_on_mesh(s4, 10, 10, envmat); 
    assert(abs(sp4.temp() - s4.temp_dist()[5].temp()) < ROUNDING);   

    InfCylinder ic4 = InfCylinder(0.1,mat,t_init);
    InfCylinderPoint icp4 = InfCylinderPoint(0.05,10);
    temp_at_point(ic4, icp4, envmat);
    assert(abs(icp4.temp()-499.973)<ROUNDING);
    temp_on_mesh(ic4, 10, 10, envmat); 
    assert(abs(icp4.temp() - ic4.temp_dist()[5].temp()) < ROUNDING);    

    PlaneWall w4 = PlaneWall(0.1, mat, t_init);
    PlaneWallPoint wp4 = PlaneWallPoint(0.05f, 10);
    temp_at_point(w4, wp4, envmat);
    assert(abs(wp4.temp()-499.973)<ROUNDING);
    temp_on_mesh(w4, 10, 10, envmat);
    assert(abs(wp4.temp() - w4.temp_dist()[5].temp()) < ROUNDING); 

    //Plot testing
    Sphere s5 = Sphere(0.0225, "egg", 293);
    vector<SpherePoint> p;
    SpherePoint sp5 = SpherePoint(0.0,0);
    SpherePoint sp6 = SpherePoint(0.0225,0);
    p.push_back(sp5);
    p.push_back(sp6);

    envmat = EnvMat("water", 358);
    plot(s5, 0, 3000, .1, envmat);
}
