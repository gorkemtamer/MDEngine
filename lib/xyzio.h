#ifndef XYZIO_H
#define XYZIO_H

#include <string>
#include "rng.h"
#include "frame.h"


//Read one frame.
//Throw std::runtime_error if file is not opening,
//malformed, contains unknown element or particles
//closer than the min distance
XYZFrame xyz_read(const std::string& path);

//write a xyz file of a given geometry
void xyz_write(const std::string& path, const XYZFrame& f);

//generate a xyz for n randomly distributed particles in a normalized box (1 1 1),
//scale coordinates and box by `scale`, write it to `path` and return the frame
//velocities drawn from Maxwell-Boltzmann distribution at T
//and center-of-mass velocity is removed
//All particles are 'Ar'
// Throws std::invalid_argument for n <= 0 or scale <= 0, and
// std::runtime_error if the box is too small to place all particles.
XYZFrame RandomFrame(int n, double scale, Boundary bc, double T, Rng& rng);
#endif