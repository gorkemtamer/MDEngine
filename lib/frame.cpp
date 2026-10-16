#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "frame.h"
#include "forcefield.h"

// Simple check for particle count (placeholder)
static int checked_count(int n){
    if (n < 0) throw std::invalid_argument("XYZFrame: negative particle count");
    return n;
}

XYZFrame::XYZFrame(int n, Vec3 box, Boundary boundary)
    :box_(box), 
    bc_(boundary),
    type_(checked_count(n), -1),
    x_(n), y_(n), z_(n),
    vx_(n), vy_(n), vz_(n),
    fx_(n), fy_(n), fz_(n)
{
    if (box.x <= 0 || box.y <= 0 || box.z <= 0)
        throw std::invalid_argument("XYZFrame: box lengths must be positive");
}

void XYZFrame::set_pos(int i, Vec3 r){
    x_[i]=r.x; y_[i]=r.y; z_[i]=r.z;
    forces_valid_=false;
}

void XYZFrame::set_element(int i, int element_index){
    type_[i] = element_index;
    forces_valid_=false;
}

void XYZFrame::clear_forces(){
    std::fill(fx_.begin(), fx_.end(), 0.0);
    std::fill(fy_.begin(), fy_.end(), 0.0);
    std::fill(fz_.begin(), fz_.end(), 0.0);
}

double XYZFrame::KineticEnergy() const {
    double ke = 0.0;
    for (int i=0; i < n(); ++i) ke += 0.5 * mass(i) *norm2(vel(i));
    return ke;
}

double XYZFrame::PotentialEnergy(const ForceField& ff){
    //forces and energy come from the same computation
    if (!forces_valid_) ff.Compute(*this);
    return epot_;
}

double XYZFrame::TotalMass() const {
    double m = 0.;
    for (int i=0; i < n(); ++i) m += mass(i);
    return m;
}

Vec3 XYZFrame::TotalMomentum() const {
    Vec3 p{};
    for (int i=0; i<n(); ++i) p += mass(i)*vel(i);
    return p;
}

Vec3 XYZFrame::CenterOfMass() const {
    Vec3 s{};
    for (int i=0; i<n(); ++i) s += mass(i) * pos(i);
    return s / TotalMass();
}

Vec3 XYZFrame::MinImage(Vec3 d) const {
    return bc_ == Boundary::Periodic ? min_image(d, box_) : d;
}

Vec3 XYZFrame::Wrap(Vec3 r) const {
    if (bc_ == Boundary::Open) return r;
    return {r.x - box_.x * std::floor(r.x / box_.x),
    r.y - box_.y * std::floor(r.y / box_.y),
    r.z - box_.z * std::floor(r.z / box_.z)};
}