#ifndef XYZIO_H
#define XYZIO_H

#include "rng.h"
#include "frame.h"

//free memory
void xyz_free(XYZFrame *f);

//write a xyz file of a given geometry
int xyz_write(const char *path, const XYZFrame *f);

//generate a xyz for n randomly distributed particles in a normalized box (1 1 1),
//scale coordinates and box by `scale`, write it to `path` and return the frame
XYZFrame* xyz_random(const char *path, int n, double scale, Rng *rng);

XYZFrame* xyz_read(const char *path);

#endif