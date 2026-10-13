#pragma once
#include <cmath>

/*Helper functions for vector operations in general
to avoid repeated function in separate sources*/

/*vector definition*/
struct Vec3{    
    double x = 0, y = 0, z = 0;

    Vec3& operator+=(const Vec3& b){x += b.x; y += b.y; z += b.z; return *this;}
    Vec3& operator-=(const Vec3& b){x -= b.x; y -= b.y; z -= b.z; return *this;}
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this;}
    Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this;}
};

/*vector operations*/
inline Vec3 operator+(Vec3 a, const Vec3& b) {return a += b;}// a+b
inline Vec3 operator-(Vec3 a, const Vec3& b) {return a -= b;}// a-b
inline Vec3 operator-(Vec3 a) {return { -a.x, -a.y, -a.z};}  //a->-a
inline Vec3 operator*(Vec3 a, double s) {return a *= s;} // a->s*a 
inline Vec3 operator*(double s, Vec3 a) {return a *= s;} // a -> a*s (commute)
inline Vec3 operator/(Vec3 a, double s) {return a /= s;} // a -> a/s

/*dot-cross products and norms*/
inline double dot(const Vec3& a, const Vec3& b) {return a.x*b.x + a.y*b.y + a.z*b.z;}
inline double norm2(const Vec3& a) {return dot(a, a);}
inline double norm(const Vec3& a) {return sqrt(norm2(a));}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

/*distance to nearest periodic image*/
inline Vec3 min_image(Vec3 d, const Vec3& L){
    d.x -= L.x * std::round(d.x / L.x);
    d.y -= L.y * std::round(d.y / L.y);
    d.z -= L.z * std::round(d.z / L.z);
    return d;
}
