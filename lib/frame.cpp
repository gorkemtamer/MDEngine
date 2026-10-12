#include <stdexcept>
#include "frame.h"

// Simple check for particle count (placeholder)
static int checked_count(int n){
    if (n < 0) throw std::invalid_argument("XYZFrame: negative particle count");
    return n;
}

XYZFrame::XYZFrame(int n_, Vec3 box_, Boundary boundary_)
    : n(checked_count(n_)), 
    box(box_), 
    boundary(boundary_),
    type(n_, -1),
    x(n_), y(n_), z(n_),
    vx(n_), vy(n_), vz(n_),
    fx(n_), fy(n_), fz(n_)
{
    if (box.x <= 0 || box.y <= 0 || box.z <= 0)
        throw std::invalid_argument("XYZFrame: box lengths must be positive");
}
