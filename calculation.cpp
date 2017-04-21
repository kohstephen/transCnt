#include <math.h>
#include <tr1/cmath>
#include "calculation.h"
#include <iostream>

extern const float PI = 3.14159265;
float RATIO = 0.003;
float PRECESION = 0.001;
vector<float> J0zeros = {0, 2.4048, 5.5201, 8.6537, 11.7915, 14.9309};

float biot(float h, float k, float x){
	return h*x/k;
}

float fourier(float alpha, float time, float x){
	return alpha*time/(x*x);
}

float calculate_alpha(float k, float density, float c){
    return k/(density*c);
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

float infinitecylinder_lumped_cap_at_time(InfiniteCylinder cylinder, float density, float h, 
    float c, float time, float t_init, float t_inf){
    float theta = exp(-h*2*time/(density*cylinder.getRadius()*c));
    return theta*(t_init-t_inf)+t_inf;
}

float planewall_lumped_cap_at_time(PlaneWall w, float density, float h, float c, float time, float t_init, float t_inf){
    float theta = exp(-h*time/(density*w.length()*c));
    return theta*(t_init-t_inf)+t_inf;
}

float sphere_lumped_cap_at_time(Sphere s, float density, float h, float c, float time, float t_init, float t_inf){
	float theta = exp(-h*3*time/(density*s.radius()*c));
	return theta*(t_init-t_inf)+t_inf;
}

float infinitecylinder_one_term_at_time_at_point(float fourier, float biot, float r, float r0, float t_init, float t_inf){
    float zeta = cylinder_solve_for_zeta(biot,1);
    float j0 = std::tr1::cyl_bessel_j(0,zeta);
    float j1 = std::tr1::cyl_bessel_j(1,zeta);
    float c1 = 2*j1/(zeta*(j0*j0+j1*j1));
    float theta = c1*exp(-zeta*zeta*fourier)*std::tr1::cyl_bessel_j(0,zeta*r/r0);
    return theta*(t_init-t_inf)+t_inf;
}

float planewall_one_term_at_time_at_point(float fourier, float biot, float x, float L, float t_init, float t_inf){
    float zeta = planewall_solve_for_zeta(biot,1);
    float c1 = 4.0f*sin(zeta)/(2.0f*zeta+sin(2.0f*zeta));
    float theta = c1*exp(-zeta*zeta*fourier)*cos(zeta*x/L);
    return theta*(t_init-t_inf)+t_inf;
}

float sphere_one_term_at_time_at_point(float fourier, float biot, float r, float r0, float t_init, float t_inf){
    float zeta = sphere_solve_for_zeta(biot,1); //can we solve this directly? (biot number)
    float c1 = 4.0f*(sin(zeta)-zeta*cos(zeta))/(2.0f*zeta-sin(2.0f*zeta));
    float temp = zeta*r/r0;
    float theta = c1*exp(-zeta*zeta*fourier)*(1/temp)*sin(temp);
    return theta*(t_init-t_inf)+t_inf;
}

float infinitecylinder_multiple_term_at_time_at_point(float fourier, float biot, float r, float r0, float t_init, float t_inf){
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
    return theta*(t_init-t_inf)+t_inf;
}

float planewall_multiple_term_at_time_at_point(float fourier, float biot, float x, float L, float t_init, float t_inf){
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
    return theta*(t_init-t_inf)+t_inf;
}

float sphere_multiple_term_at_time_at_point(float fourier, float biot, float r, float r0, float t_init, float t_inf){
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
    return theta*(t_init-t_inf)+t_inf;
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

float semi_infinite_at_time_at_point(float x, float alpha, float time, float h, float k, float t_init, float t_inf){
    float y = sqrt(alpha*time);
    float z = x/(2*y);
    float theta = erfc(z) - exp(h*x/k + h*h*alpha*time/(k*k))*erfc(z+h*y/k);
    return theta*(t_init-t_inf)+t_inf;
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

void temp_at_point(Sphere &s, SpherePoint &p, string envmat, Temp t_inf){

    float r0 = s.radius();
    string mat = s.mat();
    Temp t_init = s.t_init();
    Loc r = p.sphere_loc();
    Secs time = p.time();
    
	//heat transfer coefficient (units: W/m^2K)
    float h = get_h(envmat);
    //conduction coefficient (units: W/mK)
    float k = get_k(mat, t_init);
    float bi = biot(h, k, r0);
	//Lumped Capacitance. Use this when Bi < 0.1
    float density = get_density(mat);
    float c = get_c(mat, t_init);
    if(bi<0.1){
		p.temp( sphere_lumped_cap_at_time(s, density, h, c, time, t_init, t_inf));
		return;
	}

	//thermal diffusivity (units: m^2/s)
	float alpha = calculate_alpha(k, density, c);
	float fo = fourier(alpha, time, r0);
	//One-Term Approximation. Use this when Fo > 0.2
	if(fo > 0.2){
        p.temp( sphere_one_term_at_time_at_point(fo, bi, r, r0, t_init, t_inf));
        return;
	}
    
    //Multiple-Term Approximation.
    if(fo > 0.05){
        p.temp(sphere_multiple_term_at_time_at_point(fo, bi, r, r0, t_init, t_inf));
        return;
    }
    
    //semi-infinite approximation
    
    p.temp(semi_infinite_at_time_at_point(r0-r, alpha, time, h, k, t_init, t_inf)); 
}


void temp_at_point(PlaneWall &w, PlaneWallPoint &p, string envmat, Temp t_inf){

    float L = w.length();
    string mat = w.mat();
    Temp t_init = w.t_init();
    Loc x = p.rect_loc();
    Secs time = p.time();
    
    //heat transfer coefficient (units: W/m^2K)
    float h = get_h(envmat);
    //conduction coefficient (units: W/mK)
    float k = get_k(mat, t_init);
    float bi = biot(h, k, L);
    //Lumped Capacitance. Use this when Bi < 0.1
    float density = get_density(mat);
    float c = get_c(mat, t_init);
    if(bi<0.1){
        p.temp( planewall_lumped_cap_at_time(w, density, h, c, time, t_init, t_inf) );
        return;
    }

    //thermal diffusivity (units: m^2/s)
    float alpha = calculate_alpha(k, density, c);
    float fo = fourier(alpha, time, L);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        p.temp( planewall_one_term_at_time_at_point(fo, bi, x, L, t_init, t_inf));
        return;
    }
    
    //Multiple-Term Approximation.
    if(fo > 0.05){
        p.temp (planewall_multiple_term_at_time_at_point(fo, bi, x, L, t_init, t_inf));
        return;
    }
    
    //semi-infinite approximation
    
    p.temp( semi_infinite_at_time_at_point(L-x, alpha, time, h, k, t_init, t_inf) );
}


float temp_at_time_at_point(InfiniteCylinder cylinder, string mat, string envmat, 
    float r, float time, float t_init, float t_inf){

    float r0 = cylinder.getRadius();

    //heat transfer coefficient (units: W/m^2K)
    float h = get_h(envmat);
    //conduction coefficient (units: W/mK)
    float k = get_k(mat, t_init);
    float bi = biot(h, k, r0);
    //Lumped Capacitance. Use this when Bi < 0.1
    float density = get_density(mat);
    float c = get_c(mat, t_init);
    if(bi<0.1){
        return infinitecylinder_lumped_cap_at_time(cylinder, density, h, c, time, t_init, t_inf);
    }

    //thermal diffusivity (units: m^2/s)
    float alpha = calculate_alpha(k, density, c);
    float fo = fourier(alpha, time, r0);
    //One-Term Approximation. Use this when Fo > 0.2
    if(fo > 0.2){
        return infinitecylinder_one_term_at_time_at_point(fo, bi, r, r0, t_init, t_inf);
    }
    
    //Multiple-Term Approximation.
    if(fo > 0.05){
        return infinitecylinder_multiple_term_at_time_at_point(fo, bi, r, r0, t_init, t_inf);
    }
    
    //semi-infinite approximation
    
    return semi_infinite_at_time_at_point(r0-r, alpha, time, h, k, t_init, t_inf); 
}

/*

float temp_at_time_at_point(RectangularParallelepiped rp, string mat, string envmat, 
    float x1, float x2, float x3, float time, float t_init, float t_inf){
    PlaneWall pl1 = PlaneWall(rp.getL1());
    PlaneWall pl2 = PlaneWall(rp.getL2());
    PlaneWall pl3 = PlaneWall(rp.getL3());
    float t1 = temp_at_time_at_point(pl1, mat, envmat, x1, time, t_init, t_inf);
    float t2 = temp_at_time_at_point(pl2, mat, envmat, x2, time, t_init, t_inf);
    float t3 = temp_at_time_at_point(pl3, mat, envmat, x3, time, t_init, t_inf);
    float theta1 = (t1-t_inf)/(t_init-t_inf);
    float theta2 = (t2-t_inf)/(t_init-t_inf);
    float theta3 = (t3-t_inf)/(t_init-t_inf);
    return theta1*theta2*theta3*(t_init-t_inf)+t_inf;
}

float temp_at_time_at_point(Cylinder cylinder, string mat, string envmat, 
    float r, float x, float time, float t_init, float t_inf){
    InfiniteCylinder ic = InfiniteCylinder(cylinder.getRadius());
    PlaneWall pl = PlaneWall(cylinder.getL());
    float t1 = temp_at_time_at_point(ic, mat, envmat, r, time, t_init, t_inf);
    float t2 = temp_at_time_at_point(pl, mat, envmat, x, time, t_init, t_inf);
    float theta1 = (t1-t_inf)/(t_init-t_inf);
    float theta2 = (t2-t_inf)/(t_init-t_inf);
    return theta1*theta2*(t_init-t_inf)+t_inf;
}

float temp_at_time_at_point(InfiniteRectangularBar irb, string mat, string envmat, 
    float x1, float x2, float time, float t_init, float t_inf){
    PlaneWall pl1 = PlaneWall(irb.getL1());
    PlaneWall pl2 = PlaneWall(irb.getL2());
    float t1 = temp_at_time_at_point(pl1, mat, envmat, x1, time, t_init, t_inf);
    float t2 = temp_at_time_at_point(pl2, mat, envmat, x2, time, t_init, t_inf);
    float theta1 = (t1-t_inf)/(t_init-t_inf);
    float theta2 = (t2-t_inf)/(t_init-t_inf);
    return theta1*theta2*(t_init-t_inf)+t_inf;
}


*/