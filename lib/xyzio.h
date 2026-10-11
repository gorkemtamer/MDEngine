#ifndef XYZIO_H
#define XYZIO_H

#include "rng.h"

#define LABEL_LEN 8 //number of columns

// ---------------------------------------------------------------------
// The "dataframe": one array per column, all sharing the same length n.
// ---------------------------------------------------------------------
typedef struct {
    int n;                 // number of rows (particles)
    char comment[256];     // the XYZ comment line (line 2) specifics to the system if required
    double lattice[3];     // box lengths Lx, Ly, Lz (line 3); the box spans [0,Lx] x [0,Ly] x [0,Lz]
    int periodic;          // boundary conditions (line 3): 1 = periodic (PBC), 0 = open (OBC),
                           // -1 = not set yet (xyz_alloc); checks refuse an unset frame
    char (*label)[LABEL_LEN];  // n x LABEL_LEN char array 
    double *x, *y, *z;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;  // forces (not stored in files), filled by ComputeForces
    double epot;           // potential energy, filled by ComputeForces
    int forces_valid;      // 1 if fx, fy, fz, epot belong to the current positions;
                           // set to 0 whenever positions are changed by hand
} XYZFrame;

//memory allocation for a geometry
XYZFrame* xyz_alloc(int n);

//free memory
void xyz_free(XYZFrame *f);

//write a xyz file of a given geometry
int xyz_write(const char *path, const XYZFrame *f);

//generate a xyz for n randomly distributed particles in a normalized box (1 1 1),
//scale coordinates and box by `scale`, write it to `path` and return the frame
XYZFrame* xyz_random(const char *path, int n, double scale, Rng *rng);

XYZFrame* xyz_read(const char *path);

#endif