#ifndef CONSTANT_H
#define CONSTANT_H
#include <map>
#include <vector>
using namespace std;

/**
 * Some type alias
 */
using Kelvin = float;
using Loc = float;
using Secs = float;
using Dim = float;

/**
 * Map from material to a vector of <temparature,k> pair
 * k is the conduction coefficient (units: W/mK)
 */
static map<string, vector<pair<Kelvin,float> > > mat_to_k {
    {"al", { {100, 302}, {200, 237}, {300, 237}, {400, 240}, {600, 231}, {800, 218} } },
    {"egg", { {100, .56}, {500, .56}} },
    {"st", { {300, 60.5}, {400, 56.7}, {600, 48}, {800, 39.2}, {1000, 30} } },
    {"304", { {100, 9.2}, {200, 12.6}, {300, 14.9}, {400, 16.6}, {600, 19.8}, {800, 22.6}, {1000, 25.4},
    {1200, 28}, {1500, 31.7}} } };

/**
 * Map from material to a vector of <temparature, c> pair
 * c is 'the ratio of the heat added to (or removed from) an object
 * to the resulting temperature change” (Wikipedia)  (units: J/K)
 */
static map<string, vector<pair<Kelvin,float> > > mat_to_c {
    {"al", { {100, 482}, {200, 798}, {300, 903}, {400, 949}, {600, 1033}, {800, 1146} } },
    {"egg", { {100, 3400}, {500, 3400}} },
    {"st", { {300, 434}, {400, 487}, {600, 559}, {800, 685}, {1000, 1169} } }, 
    {"304", { {100, 272}, {200, 402}, {300, 477}, {400, 515}, {600, 557}, {800, 582}, {1000, 611},
    {1200, 640}, {1500, 682}} }
};

/**
 * Map from material to h
 * h is the heat transfer coefficient (units: W/m2K)
 */
static map<string, float> mat_to_h {
    {"air", 100},
    {"air2", 10},
    {"water", 1000}, 
    {"oil", 500}
};

/**
 * Map from material to p
 * p is the density of the object (units: kg/m3)
 */
static map<string, float> mat_to_p {
    {"al", 2702},
    {"egg", 1202.233},
    {"st", 7854}, 
    {"304", 7900}
};


/**
 * Calculate k from mat and temparature.
 * We assume k is linear between two temparatures.
 */ 
static float get_k(string mat, Kelvin temp){
    // left to do: check if the requested temp is in range
    vector<pair<float,float> > pairs = mat_to_k[mat];
    pair<float,float> prev;
    for(auto it = pairs.begin(); it != pairs.end(); ++it){
        if( it->first < temp ){
            prev = *it;
        }else{
            return prev.second+(it->second - prev.second) * (temp-prev.first)/(it->first-prev.first);
        }
    }
    return -1;
}

/**
 * Calculate c from mat and temparature.
 * We assume c is linear between two temparatures.
 */ 
static float get_c(string mat, float temp){
    vector<pair<float,float>> pairs = mat_to_c[mat];
    pair<float,float> prev;
    for(auto it = pairs.begin(); it != pairs.end(); ++it){
        if( it->first < temp ){
            prev = *it;
        }else{
            return prev.second+(it->second - prev.second) * (temp-prev.first)/(it->first-prev.first);
        }
    }
    return -1;
}

/**
 * Get corresponding h of the given envmat.
 */ 
static float get_h(string envmat){
    return mat_to_h[envmat];
}

/**
 * Get correspoinding p of the given mat.
 */
static float get_p(string mat){
    return mat_to_p[mat];
}

#endif
