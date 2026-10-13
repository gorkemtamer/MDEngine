#pragma once
#include <string>
#include <vector>
#include "vec3.h"
#include "boundary.h"
#include "elements.h"

inline constexpr int LABEL_LEN = 8; //number of columns in file + 1

//--------------------------------------------------
/*
One snapshot of the system. All per-particle arrays
have length n.
*/
//--------------------------------------------------

class XYZFrame {
public:
    XYZFrame(int n, Vec3 box, Boundary bc); 

    //#particles, box, boundary conditions
    int n() const {return static_cast<int>(type_.size());}
    Vec3 box() const {return box_;}
    Boundary boundary() const {return bc_;}

    // per particle
    Vec3 pos(int i) const {return {x_[i],y_[i],z_[i]};}
    Vec3 vel(int i) const {return {vx_[i], vy_[i], vz_[i]};}
    Vec3 force(int i) const {return {fx_[i], fy_[i], fz_[i]};}
    void set_pos(int i, Vec3 r);
    void set_vel(int i, Vec3 v) {vx_[i]=v.x; vy_[i]=v.y;vz_[i]=v.z;}
    
    int type(int i) const {return type_[i];} //element index
    const Element& element(int i) const { return ::element(type_[i]);}
    const char* label(int i) const {return element(i).label;}
    double mass(int i) const {return element(i).mass; }
    void set_element(int i, int element_index);
    
    //force related part
    void clear_forces();
    void add_force(int i, Vec3 F) {fx_[i] += F.x; fy_[i] += F.y; fz_[i] += F.z;}
    void set_potential_energy(double U) {epot_= U; forces_valid_ = true;}
    double PotentialEnergy() const {return epot_;}
    bool forces_valid() const {return forces_valid_;}
    
    std::string comment; //comment in the input file (line 2)

    //rest of the physics
    double KineticEnergy() const;
    double TotalMass() const;
    Vec3 TotalMomentum() const;
    Vec3 CenterOfMass() const; 
    Vec3 MinImage(Vec3 d) const; //plain d for OBC, nearest image for PBC
    Vec3 Wrap(Vec3 r) const; //into the box for PBC, nothing for OBC

private:
    Vec3 box_;
    Boundary bc_;
    std::vector<int> type_; // element index list 
    std::vector<double> x_,y_,z_,vx_,vy_,vz_,fx_,fy_,fz_;
    double epot_ = 0.0;
    bool forces_valid_ = false; //true if fx,fy,fz,epot belong to the current positions
};