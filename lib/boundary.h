#pragma once

//---------------------------------------------------
/*
Boundary condition definition:
Open -> no periodicity, direct computation of the distances
Periodic -> the box is periodic, a particle also appears
        on the other side if the border is crossed,
        distances use the nearest periodic image (min_image in vec3.h)
Defined in .xyz files as the flag periodic 0 (OBC) or periodic 1 (PBC)
*/ 
//---------------------------------------------------

enum class Boundary {Open, Periodic};

//value from .xyz file
constexpr int boundary_flag(Boundary b){
    return b == Boundary::Periodic ? 1 : 0;
}

//flag read from the file
constexpr bool boundary_from_flag(int flag, Boundary &out){
    if (flag==0) {out=Boundary::Open; return true;}
    if (flag==1) {out=Boundary::Periodic; return true;}
    return false;
}

constexpr const char* boundary_name(Boundary b){
    return b == Boundary::Periodic ? "periodic" : "open";
}