#include <algorithm>
#include <stdexcept>
#include "frame.h"

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