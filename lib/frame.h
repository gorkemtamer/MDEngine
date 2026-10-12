#pragma once
#include <string>
#include <vector>
#include "vec3.h"
#include "boundary.h"
#include "elements.h"

#define LABEL_LEN 8 

//--------------------------------------------------
/*
One snapshot of the system. All per-particle arrays
have length n.
*/
//--------------------------------------------------

class XYZFrame{
public:

    XYZFrame(int n, Vec3 box, Boundary boundary);

    int n; //#particles
    std::string comment; //comment line in xyz file
    Vec3 box; //box sizes
    Boundary boundary; //periodic 0:OBC, periodic 1:PBC

    std::vector<int> type; //element index list
    std::vector<double> x,y,z;
    std::vector<double> vx,vy,vz;
    std::vector<double> fx,fy,fz; 

    double epot = 0.0;
    int forces_valid = 0; // 1 if fx,fy,fz,epot belong to the current position

    //element symbol of particle i
    const char* label(int i) const {return element(type[i]).label;}

};

//redefinition of XYZFrame as a class
//compared to struct XYZFrame in C



// class XYZFrame {
// public:
//     XYZFrame(int n, Vec3 box, Boundary bc); 

//     //#particles, box, boundary conditions
//     int n() const {return static_cast<int>(type_.size());}
//     Vec3 box() const {return box_;}
//     Boundary boundary() const {return bc_;}

//     // per particle
//     Vec3 pos(int i) const {return {x_[i],y_[i],z_[i]};}
//     Vec3 vel(int i) const {return {vx_[i], vy_[i], vz_[i]};}
//     Vec3 force(int i) const {return {fx_[i], fy_[i], fz_[i]};}
//     void set_pos(int i, Vec3 r);
//     void set_vel(int i, Vec3 v);
//     void set_element(int i, int element_index);
//     const Element& element(int i) const; //-> label, sigma, eps, mass
//     double mass(int i) const;

//     //physics part
//     double KineticEnergy() const;
//     Vec3 TotalMomentum() const;
//     Vec3 CenterOfMass() const;
//     Vec3 MinImage(Vec3 d) const; //plane d for OBC, nearest image for PBC
//     Vec3 wrap(Vec3 r) const; //into the box for PBC, nothing for OBC

//     //forces and potential energy
//     void clear_forces();
//     void add_force(int i, Vec3 F);
//     void set_potential_energy(double U); //also marks the forces as current
//     double PotentialEnergy() const;
//     bool forces_valid() const;

//     std::string comment; //comment in the input file
// private:
//     Vec3 box_;
//     Boundary bc_;
//     std::vector<int> type_; // element index list 
//     std::vector<double> x_,y_,z_,vx_,vy_,vz_,fx_,fy_,fz_;
//     double epot_ = 0.0;
//     bool forces_valid_ = false;
// };