#include "constant.h"
// a Point is a concept that holds a time (built-in duration type). 

template<typename T>
concept bool Point = requires(T p, Temp t) {
        { p.temp() } -> Temp;
        { p.time() } -> Secs;
        { p.temp(t) };
};


// There will be 6 concepts refined from the Point concept:

template<typename T>
concept bool PlaneWallPoint = Point<T> && requires(T p) {
        { p.rect_loc() } -> Loc;
};

template<typename T>
concept bool InfCylPoint = Point<T> && requires(T p) {
        { p.cyl_loc() } -> Loc;
};

template<typename T>
concept bool SpherePoint = Point<T> && requires(T p) {
        { p.sphere_loc() } -> Loc;
};

template<typename T>
concept bool InfRectBarPoint = PlaneWallPoint<T> && requires(T p) {
        { p.rect_loc2() } -> Loc;
};

template<typename T>
concept bool RectBarPoint = InfRectBarPoint<T> && requires(T p) {
        { p.rect_loc3() } -> Loc;
};

template<typename T>
concept bool ShortCylPoint = InfCylPoint<T> && requires(T p) {
        { p.rect_loc() } -> Loc;
};