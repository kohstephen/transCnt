#include <math.h>
#include <tr1/cmath>
#include <iostream>
#include "calculation.h"
#include <valarray>

extern const float PI = 3.14159265;
//TODO: change "RATIO" to 0.01
float RATIO = 0.003;
float PRECESION = 0.001;
float EPSILON = 0.01;
vector<float> J0zeros = {0, 2.4048, 5.5201, 8.6537, 11.7915, 14.9309};

float kelvin_to_fahrenheit(Kelvin k){
    return k*9.0/5.0 - 459.67;
}

float kelvin_to_celcius(Kelvin k){
    return k - 273.15;
}

float biot(float h, float k, float x){
    return h*x/k;
}

float fourier(float alpha, float time, float x){
    return alpha*time/(x*x);
}


float theta_to_temp(float theta, float t_init, float t_inf){
    return theta*(t_init-t_inf)+t_inf;
}

// only support n <= 6
float cylinder_solve_for_zeta(float biot, int n){
    float min = J0zeros[n-1];
    float max = J0zeros[n];
    float x;
    float fx;
    while(true){
        x = (min+max)/2;
        fx = x*std::tr1::cyl_bessel_j(1,x)/std::tr1::cyl_bessel_j(0,x);
        if(abs(fx-biot)<PRECESION) return x;
        if(fx<biot) min = x;
        else max = x;
    }
    return -1;
}

float planewall_solve_for_zeta(float biot, int n){
    // [0, pi/2] [pi*(n-3/2), pi*(n-1/2)]
    //ζtan(ζ) = Bi
    float min = n==1? 0:(n-1.5)*PI;
    float max = (n-0.5)*PI;
    float x;
    float fx;
    while(true){
        x = (min+max)/2;
        fx = x*tan(x);
        if(abs(fx-biot)<PRECESION) return x;
        if(fx<biot) min = x;
        else max = x;
    }
    //return 1.142;
    return -1;
}

//newton-raphson
float sphere_solve_for_zeta(float biot, int n){
    // [(n-1)*pi, n*pi]
    //ζcot(ζ) = 1-Bi
    float rh = 1-biot;
    float min = (n-1)*PI;
    float max = n*PI;
    float x;
    float fx;
    while(true){
        x = (min+max)/2;
        fx = x/tan(x);
        if(abs(fx-rh)<PRECESION) return x;
        if(fx<rh) max = x;
        else min = x;
    }
    //return 1.142;
    return -1;
}

float infinitecylinder_lumped_cap_at_time(float r0, float density, float h, float c, float time){
    float theta = exp(-h*2*time/(density*r0*c));
    return theta;
}

float planewall_lumped_cap_at_time(float L, float density, float h, float c, float time){
    float theta = exp(-h*time/(density*L*c));
    return theta;
}

float sphere_lumped_cap_at_time(float r0, float density, float h, float c, float time){
    float theta = exp(-h*3*time/(density*r0*c));
    return theta;
}

float infinitecylinder_one_term_at_time_at_point(float fourier, float biot, float r, float r0){
    float zeta = cylinder_solve_for_zeta(biot,1);
    float j0 = std::tr1::cyl_bessel_j(0,zeta);
    float j1 = std::tr1::cyl_bessel_j(1,zeta);
    float c1 = 2*j1/(zeta*(j0*j0+j1*j1));
    float theta = c1*exp(-zeta*zeta*fourier)*std::tr1::cyl_bessel_j(0,zeta*r/r0);
    return theta;
}

float planewall_one_term_at_time_at_point(float fourier, float biot, float x, float L){
    float zeta = planewall_solve_for_zeta(biot,1);
    float c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
    float theta = c1*exp(-zeta*zeta*fourier)*cos(zeta*x/L);
    return theta;
}

float sphere_one_term_at_time_at_point(float fourier, float biot, float r, float r0){
    float zeta = sphere_solve_for_zeta(biot,1); //can we solve this directly? (biot number)
    float c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
    float temp = zeta*r/r0;
    float theta = c1*exp(-zeta*zeta*fourier)*(1/temp)*sin(temp);
    return theta;
}

float infinitecylinder_multiple_term_at_time_at_point(float fourier, float biot, float r, float r0){
    float theta = 0;
    int n = 1;
    float zeta_n, zeta_n_1, j0, j1, c_n, c_n_1;
    while(true){
        if(n>1){
            zeta_n_1 = zeta_n;
            c_n_1 = c_n;
        }
        zeta_n = cylinder_solve_for_zeta(biot, n);
        j0 = std::tr1::cyl_bessel_j(0,zeta_n);
        j1 = std::tr1::cyl_bessel_j(1,zeta_n);
        c_n = 2*j1/(zeta_n*(j0*j0+j1*j1));
        if(n>1){
            float ratio = abs(c_n*exp(-zeta_n*zeta_n*fourier)/(c_n_1*exp(-zeta_n_1*zeta_n_1*fourier)));
            if(ratio < RATIO) break;
        }
        theta += c_n*exp(-zeta_n*zeta_n*fourier)*std::tr1::cyl_bessel_j(0, zeta_n*r/r0);
        n += 1;
    }
    return theta;
}

float planewall_multiple_term_at_time_at_point(float fourier, float biot, float x, float L){
    float theta = 0;
    int n = 1;
    float zeta_n, zeta_n_1, c_n, c_n_1;
    while(true){
        if(n>1){
            zeta_n_1 = zeta_n;
            c_n_1 = c_n;
        }
        zeta_n = planewall_solve_for_zeta(biot, n);
        c_n = 4.0f*sin(zeta_n)/(2.0f*zeta_n+sin(2.0f*zeta_n));
        if(n>1){
            float ratio = abs(c_n*exp(-zeta_n*zeta_n*fourier)/(c_n_1*exp(-zeta_n_1*zeta_n_1*fourier)));
            if(ratio < RATIO) break;
        }
        theta += c_n*exp(-zeta_n*zeta_n*fourier)*cos(zeta_n*x/L);
        n += 1;
    }
    return theta;
}

float sphere_multiple_term_at_time_at_point(float fourier, float biot, float r, float r0){
    // std::cout << biot << std::endl;
    float theta = 0;
    int n = 1;
    float zeta_n, zeta_n_1, c_n, c_n_1;
    while(true){
        if(n>1){
            zeta_n_1 = zeta_n;
            c_n_1 = c_n;
        }
        zeta_n = sphere_solve_for_zeta(biot, n);
        // std::cout << zeta_n << std::endl;
        c_n = 4.0f*(sin(zeta_n)-zeta_n*cos(zeta_n))/(2.0f*zeta_n-sin(2.0f*zeta_n));
        // std::cout << c_n << std::endl;
        if(n>1){
            float ratio = abs(c_n*exp(-zeta_n*zeta_n*fourier)/(c_n_1*exp(-zeta_n_1*zeta_n_1*fourier)));
            if(ratio < RATIO) break;
        }
        float temp = zeta_n*r/r0;
        theta += c_n*exp(-zeta_n*zeta_n*fourier)*(1/temp)*sin(temp);
        n += 1;
    }
    return theta;
}

/**
  void sphere_one_term_at_time(vector<float>* ret, float fourier, float biot, vector<float>& points, float r0, float t_init, float t_inf){
  float zeta; //can we solve this directly?
  float c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
  float y = c1*exp(-zeta*zeta*fourier);
  float z = zeta/r0;
  float diff = t_init-t_inf;
  int i = 0;
  for(auto it = points.begin(); it != points.end(); ++it, ++i){
  float q = z*(*it);
  float theta = y*(1/q)*sin(q);
  (*ret)[i] = theta*diff+t_inf;
  }
  }
 **/

float semi_infinite_at_time_at_point(float x, float alpha, float time, float h, float k){
    float y = sqrt(alpha*time);
    float z = x/(2*y);
    float theta = erfc(z) - exp(h*x/k + h*h*alpha*time/(k*k))*erfc(z+h*y/k);
    return 1 - theta;
}

/**
  vector<float> semi_infinite_at_time(vector<float>* ret, vector<float>& points, float alpha, float time, float h, float k, float t_init, float t_inf){
  float y = sqrt(alpha*time);
  float z = h*h*alpha*time/(k*k);
  float diff = t_init-t_inf;
  int i = 0;
  for(auto it = points.begin(); it != points.end(); ++it, ++i){
  float x = *it;
  float q = x/(2*y);
  float theta = erfc(q) - exp(h*x/k + z)*erfc(q+h*y/k);
  (*ret)[i] = theta*diff+t_inf;
  }
  }
 **/

/**
  void temp_at_time(vector<float>* ret, Sphere s, string mat, string envmat, vector<float>& points, float time, float t_init, float t_inf){
  float r0 = s.getRadius();

//heat transfer coefficient (units: W/m^2K)
float h = get_h(envmat);
//conduction coefficient (units: W/mK)
float k = get_k(mat, t_init);
float bi = biot(h, k, r0);
//Lumped Capacitance. Use this when Bi < 0.1
// saved data
float density = get_density(mat);
float c = get_c(mat, t_init);
if(bi<0.1){
float temp =  sphere_lumped_cap_at_time(s, density, h, c, time, t_init, t_inf);
for(int i = 0; i < points.size(); i++){
(*ret)[i] = temp;
}
return;
}

//thermal diffusivity (units: m^2/s)
float alpha = calculate_alpha(k, density, c);
float fo = fourier(alpha, time, r0);
//One-Term Approximation. Use this when Fo > 0.2
if(fo > 0.2){
sphere_one_term_at_time(ret, fo, bi, points, r0, t_init, t_inf);
return;
}

//Saved data.
if(fo > 0.05){

}

//semi-infinite approximation
semi_infinite_at_time(ret, points, alpha, time, h, k, t_init, t_inf);

}
 **/
 
float theta_at_point(Sphere &s, SpherePoint &p, float h){
    float r0 = s.radius();
    Loc r = p.sphere_loc();
    Secs time = p.time();
    float k = s.k(); //conduction coefficient (units: W/mK)
    float bi = biot(h, k, r0);
    float density = s.p();
    float c = s.c();

    //Lumped Capacitance. Use this when Bi < 0.1
    if(bi<0.1){
        return sphere_lumped_cap_at_time(r0, density, h, c, time);
    }

    //thermal diffusivity (units: m^2/s)
    float alpha = s.a();
    float fo = fourier(alpha, time, r0);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        return sphere_one_term_at_time_at_point(fo, bi, r, r0);
    }

    //Multiple-Term Approximation.
    if(fo > 0.05){
        return sphere_multiple_term_at_time_at_point(fo, bi, r, r0);
    }

    //semi-infinite approximation

    return semi_infinite_at_time_at_point(r0-r, alpha, time, h, k); 
}

float theta_at_point(PlaneWall &w, PlaneWallPoint &p, float h){
    float L = w.length();
    Loc x = p.rect_loc();
    Secs time = p.time();
    float k = w.k();
    float bi = biot(h, k, L);
    float density = w.p();
    float c = w.c();

    if(bi<0.1){
        return planewall_lumped_cap_at_time(L, density, h, c, time);
    }

    float alpha = w.a();
    float fo = fourier(alpha, time, L);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        return planewall_one_term_at_time_at_point(fo, bi, x, L);
    }

    //Multiple-Term Approximation.
    if(fo > 0.05){
        return planewall_multiple_term_at_time_at_point(fo, bi, x, L);
    }

    //semi-infinite approximation
    return semi_infinite_at_time_at_point(L-x, alpha, time, h, k);
}

float theta_at_point(InfCylinder &icyl, InfCylinderPoint &p, float h){
    float r0 = icyl.radius();
    Loc r = p.cyl_loc();
    Secs time = p.time();
    float k = icyl.k();
    float bi = biot(h, k, r0);
    float density = icyl.p();
    float c = icyl.c();
    if(bi<0.1){
        return infinitecylinder_lumped_cap_at_time(r0, density, h, c, time);
    }

    //thermal diffusivity (units: m^2/s)
    float alpha = icyl.a();
    float fo = fourier(alpha, time, r0);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        return infinitecylinder_one_term_at_time_at_point(fo, bi, r, r0);
    }

    //Multiple-Term Approximation.
    if(fo > 0.05){
        return infinitecylinder_multiple_term_at_time_at_point(fo, bi, r, r0);
    }

    //semi-infinite approximation

    return semi_infinite_at_time_at_point(r0-r, alpha, time, h, k);
}

void temp_at_point(Sphere &s, SpherePoint &p, EnvMat &envmat){
    p.temp(theta_to_temp(theta_at_point(s, p, envmat.h()), s.t_init(), envmat.t_inf()));
}

void temp_at_point(PlaneWall &w, PlaneWallPoint &p, EnvMat &envmat){
    p.temp(theta_to_temp(theta_at_point(w, p, envmat.h()), w.t_init(), envmat.t_inf()));
}

void temp_at_point(InfCylinder &icyl, InfCylinderPoint &p, EnvMat &envmat){
    p.temp(theta_to_temp(theta_at_point(icyl, p, envmat.h()), icyl.t_init(), envmat.t_inf()));
}

void temp_at_point(RectBar &rb, RectBarPoint &p, EnvMat &envmat){
    float k = rb.k();
    float c = rb.c();
    float density = rb.p();
    Kelvin t_init = rb.t_init();
    float h = envmat.h();
    Secs time = p.time();

    PlaneWall pl1 = PlaneWall(rb.l1(), k,c,density, t_init);
    PlaneWall pl2 = PlaneWall(rb.l2(), k,c,density, t_init);
    PlaneWall pl3 = PlaneWall(rb.l3(), k,c,density, t_init);

    PlaneWallPoint p1 = PlaneWallPoint(p.rect_loc1(),time);
    PlaneWallPoint p2 = PlaneWallPoint(p.rect_loc2(),time);
    PlaneWallPoint p3 = PlaneWallPoint(p.rect_loc3(),time);

    float theta1 = theta_at_point(pl1, p1, h);
    float theta2 = theta_at_point(pl2, p2, h);
    float theta3 = theta_at_point(pl3, p3, h);

    p.temp(theta_to_temp(theta1*theta2*theta3, t_init, envmat.t_inf()));
}


void temp_at_point(Cylinder &cyl, CylinderPoint &p, EnvMat &envmat){
    float k = cyl.k();
    float c = cyl.c();
    float density = cyl.p();
    Kelvin t_init = cyl.t_init();
    float h = envmat.h();
    Secs time = p.time();

    InfCylinder icyl = InfCylinder(cyl.radius(),k,c,density,t_init);
    PlaneWall w = PlaneWall(cyl.length(),k,c,density,t_init);

    InfCylinderPoint cylp = InfCylinderPoint(p.cyl_loc(),time);
    PlaneWallPoint wp = PlaneWallPoint(p.rect_loc(), time);

    float theta1 = theta_at_point(icyl, cylp, h);
    float theta2 = theta_at_point(w, wp, h);

    p.temp(theta_to_temp(theta1*theta2, t_init, envmat.t_inf()));
}

void temp_at_point(InfRectBar &irb, InfRectBarPoint &p, EnvMat &envmat){
    float k = irb.k();
    float c = irb.c();
    float density = irb.p();
    Kelvin t_init = irb.t_init();
    float h = envmat.h();
    Secs time = p.time();

    PlaneWall pl1 = PlaneWall(irb.l1(), k,c,density, t_init);
    PlaneWall pl2 = PlaneWall(irb.l2(), k,c,density, t_init);

    PlaneWallPoint p1 = PlaneWallPoint(p.rect_loc1(),time);
    PlaneWallPoint p2 = PlaneWallPoint(p.rect_loc2(),time);

    float theta1 = theta_at_point(pl1, p1, h);
    float theta2 = theta_at_point(pl2, p2, h);

    p.temp(theta_to_temp(theta1*theta2, t_init, envmat.t_inf()));
}

valarray<float> lumped_cap_on_mesh(PlaneWall &w, EnvMat &envmat, valarray<Loc> locs, Secs secs) {
    valarray<float> theta (exp(-1 * envmat.h()*secs/(w.p()*w.length()*w.c())), locs.size());
    return theta;
}

valarray<float> one_term_at_time_on_mesh(float fo, float bi, valarray<float> locs, Dim L) {
    float zeta = planewall_solve_for_zeta(bi, 1);
    float c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
    valarray<float> theta = c1*exp(-zeta*zeta*fo)*cos(zeta*locs/L);      
    return theta;
}

valarray<float> multiple_term_at_time_on_mesh(float fo, float bi, valarray<float> locs, Dim L) {
    valarray<float> theta (locs.size());
    int n = 1;
    float zeta_n, zeta_n_1, c_n, c_n_1;
    while(true){
        if(n>1){
            zeta_n_1 = zeta_n;
            c_n_1 = c_n;
        }
        zeta_n = planewall_solve_for_zeta(bi, n);
        c_n = 4.0f*sin(zeta_n)/(2.0f*zeta_n+sin(2.0f*zeta_n));
        if(n>1){
            float ratio = abs(c_n*exp(-zeta_n*zeta_n*fo)/(c_n_1*exp(-zeta_n_1*zeta_n_1*fo)));
            if(ratio < RATIO) break;
        }
        theta += c_n*exp(-zeta_n*zeta_n*fo)*cos(zeta_n*locs/L);
        n += 1;
    }
    return theta;     
}

valarray<float> semi_infinite_at_time_on_mesh(valarray<float> x, float alpha, float time, float h, float k){
    float y = sqrt(alpha*time);
    valarray<float> z (x/(2*y));
    // not allowed to do erfc on valarrays -- have to iterate over each element in valarray
    valarray<float> theta (x.size());
    for (int i = 0; i < theta.size(); i++) { 
        theta[i] = erfc(z[i]) - exp(h*x[i]/k + h*h*alpha*time/(k*k))*erfc(z[i]+h*y/k);
    }
    valarray<float> ones (1, x.size());
    return ones - theta;
}

valarray<float> theta_on_mesh(PlaneWall &w, Secs secs, int mesh_density, EnvMat &envmat) {
    // temp_dist(w, num_points, secs);
    int num_points = mesh_density + 1; 
    float L = w.length();
    float k = w.k();
    float h = envmat.h();
    float bi = biot(h, k, L);
    valarray<Loc> locs (num_points);  
    float incr = w.length() / mesh_density; 
    // vector<PlaneWallPoint> temp_dist;
    for (int i = 0; i < num_points; i++) {
        locs[i] = i*incr;	
    }
    if(bi<0.1){
        cout << "LUMPED" << endl;
        return lumped_cap_on_mesh(w, envmat, locs, secs);
    }

    float alpha = w.a();
    float fo = fourier(alpha, secs, L);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        cout << "ONE-TERM" << endl;
        return one_term_at_time_on_mesh(fo, bi, locs, w.length());
    }

    //Multiple-Term Approximation.
    if(fo > 0.05){
        cout << "MULTI-TERM" << endl;
        return multiple_term_at_time_on_mesh(fo, bi, locs, w.length());
    }

    //semi-infinite approximation
    cout << "SEMI-INF" << endl;
    valarray<Dim> Ls (L, locs.size());
    valarray<float> x = Ls - locs;
    return semi_infinite_at_time_on_mesh(x, alpha, secs, h, k);
}

valarray<Kelvin> theta_to_temp(valarray<float> theta, float t_init, float t_inf){
    return theta*(t_init-t_inf)+t_inf;
}

void temp_on_mesh(PlaneWall &w, Secs secs, int mesh_density, EnvMat &envmat){
    valarray<Kelvin> temps = theta_to_temp(theta_on_mesh(w, secs, mesh_density, envmat), w.t_init(), envmat.t_inf());  
    for (Kelvin i : temps) {
        cout << i << ' ';
    }
}

float avg_temp_at_time(Sphere &s, Secs time, EnvMat &envmat){
    float r0 = s.radius();
    float k = s.k();
    float h = envmat.h();
    float bi = biot(h, k, r0);
    float density = s.p();
    float c = s.c();
    float zeta, c1, theta;
    float term_cur, term_prev;

    if(bi<0.1){
        return sphere_lumped_cap_at_time(r0, density, h, c, time);
    }

    float alpha = s.a();
    float fo = fourier(alpha, time, r0);
    if(fo > 0.2){
        zeta = sphere_solve_for_zeta(bi,1);
        c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
        theta = c1*exp(-zeta*zeta*fo);
        return theta;
    }

    if(fo > 0.05){
        int n = 1;
        theta = 0;
        zeta = sphere_solve_for_zeta(bi,n);
        c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
        theta = c1*exp(-zeta*zeta*fo);
        term_prev = theta;
        ++n;
        while(true) {
            zeta = sphere_solve_for_zeta(bi,n);
            c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
            term_cur = c1*exp(-zeta*zeta*fo);
            if(term_cur/term_prev > EPSILON) {
                theta += term_cur;
                term_prev = term_cur;
                ++n;
            } else {
                return theta;
            }
        }
    }

    return s.t_init();
}


float avg_temp_at_time(PlaneWall &w, Secs time, EnvMat &envmat){
    float L = w.length();
    float k = w.k();
    float h = envmat.h();
    float bi = biot(h, k, L);
    float density = w.p();
    float c = w.c();
    float zeta, c1, theta;
    float term_cur, term_prev;

    if(bi<0.1){
        return planewall_lumped_cap_at_time(L, density, h, c, time);
    }

    float alpha = w.a();
    float fo = fourier(alpha, time, L);
    if(fo > 0.2){
        zeta = planewall_solve_for_zeta(bi,1);
        c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
        theta = c1*exp(-zeta*zeta*fo);
        return theta;
    }

    if(fo > 0.05){
        int n = 1;
        theta = 0;
        zeta = planewall_solve_for_zeta(bi,n);
        c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
        theta = c1*exp(-zeta*zeta*fo);
        term_prev = theta;
        ++n;
        while(true) {
            zeta = planewall_solve_for_zeta(bi,n);
            c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
            term_cur = c1*exp(-zeta*zeta*fo);
            if(term_cur/term_prev > EPSILON) {
                theta += term_cur;
                term_prev = term_cur;
                ++n;
            } else {
                return theta;
            }
        }
    }

    return w.t_init();
}

float avg_temp_at_time(InfCylinder &icyl, Secs time, EnvMat &envmat){
    float r0 = icyl.radius();
    float h = envmat.h();
    float k = icyl.k();
    float bi = biot(h, k, r0);
    float density = icyl.p();
    float c = icyl.c();
    float zeta, c1, theta;
    float term_cur, term_prev;
    float j0, j1;

    if(bi<0.1){
        return infinitecylinder_lumped_cap_at_time(r0, density, h, c, time);
    }

    float alpha = icyl.a();
    float fo = fourier(alpha, time, r0);
    if(fo > 0.2){
        zeta = cylinder_solve_for_zeta(bi,1);
        j0 = std::tr1::cyl_bessel_j(0,zeta);
        j1 = std::tr1::cyl_bessel_j(1,zeta);
        c1 = 2*j1/(zeta*(j0*j0+j1*j1));
        theta = c1*exp(-zeta*zeta*fo);
        return theta;
    }

    if(fo > 0.05){
        int n = 1;
        theta = 0;
        zeta = cylinder_solve_for_zeta(bi,n);
        j0 = std::tr1::cyl_bessel_j(0,zeta);
        j1 = std::tr1::cyl_bessel_j(1,zeta);
        c1 = 2*j1/(zeta*(j0*j0+j1*j1));
        theta = c1*exp(-zeta*zeta*fo);
        term_prev = theta;
        ++n;
        while(true) {
            zeta = cylinder_solve_for_zeta(bi,n);
            j0 = std::tr1::cyl_bessel_j(0,zeta);
            j1 = std::tr1::cyl_bessel_j(1,zeta);
            c1 = 2*j1/(zeta*(j0*j0+j1*j1));
            term_cur = c1*exp(-zeta*zeta*fo);
            if(term_cur/term_prev > EPSILON) {
                theta += term_cur;
                term_prev = term_cur;
                ++n;
            } else {
                return theta;
            }
        }
    }

    return icyl.t_init();
}

float avg_temp_at_time(RectBar &rb, Secs time, EnvMat &envmat){
    float k = rb.k();
    float c = rb.c();
    float density = rb.p();
    Kelvin t_init = rb.t_init();
    float h = envmat.h();

    PlaneWall pl1 = PlaneWall(rb.l1(), k,c,density, t_init);
    PlaneWall pl2 = PlaneWall(rb.l2(), k,c,density, t_init);
    PlaneWall pl3 = PlaneWall(rb.l3(), k,c,density, t_init);

    float theta1 = avg_temp_at_time(pl1, time, envmat);
    float theta2 = avg_temp_at_time(pl2, time, envmat);
    float theta3 = avg_temp_at_time(pl3, time, envmat);

    return theta1*theta2*theta3;
}

float avg_temp_at_time(Cylinder &cyl, Secs time, EnvMat &envmat){
    float k = cyl.k();
    float c = cyl.c();
    float density = cyl.p();
    Kelvin t_init = cyl.t_init();
    float h = envmat.h();

    InfCylinder icyl = InfCylinder(cyl.radius(),k,c,density,t_init);
    PlaneWall w = PlaneWall(cyl.length(),k,c,density,t_init);

    float theta1 = avg_temp_at_time(icyl, time, envmat);
    float theta2 = avg_temp_at_time(w, time, envmat);

    return theta1*theta2;
}

float avg_temp_at_time(InfRectBar &irb, Secs time, EnvMat &envmat){
    float k = irb.k();
    float c = irb.c();
    float density = irb.p();
    Kelvin t_init = irb.t_init();
    float h = envmat.h();

    PlaneWall pl1 = PlaneWall(irb.l1(), k,c,density, t_init);
    PlaneWall pl2 = PlaneWall(irb.l2(), k,c,density, t_init);

    float theta1 = avg_temp_at_time(pl1, time, envmat);
    float theta2 = avg_temp_at_time(pl2, time, envmat);

    return theta1*theta2;
}
