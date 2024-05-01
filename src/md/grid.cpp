// Author: Moses KJ Chung
// Year:   2024

#include "alfatom.h"
#include "grid.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  grid  --  sort atoms geometrically using 3D grid  //
//                                                    //
////////////////////////////////////////////////////////

// "grid" consists of a set of procedures for sorting point
// geometrically using a 3D grid for AlphaMol2

void splitGrid(AlfAtom *alfatoms, int size, real xmin, real xmax, real ymin, real ymax, real zmin, real zmax, int ncube, std::vector<int>& Nval)
{
    int Nx,Ny,Nz;
    int idx,Ix,Iy,Iz;
    real hx,hy,hz;
    real Xlim[2],Ylim[2],Zlim[2];
    real Point[3];
    real offset = 0.1;

    Xlim[0] = xmin - offset;
    Xlim[1] = xmax + offset;
    Ylim[0] = ymin - offset;
    Ylim[1] = ymax + offset;
    Zlim[0] = zmin - offset;
    Zlim[1] = zmax + offset;

    Nx = 1; Ny = 1; Nz = 1;
    if(ncube == 1) {
        Nx = 1; Ny = 1; Nz = 1;
    }
    else if(ncube == 2) {
        Nx = 2; Ny = 1; Nz = 1;
    }
    else if(ncube == 4) {
        Nx = 2; Ny = 2; Nz = 1;
    }
    else if(ncube == 8) {
        Nx = 2; Ny = 2; Nz = 2;
    }
    else if(ncube == 16) {
        Nx = 4; Ny = 2; Nz = 2;
    }
    else if(ncube == 32) {
        Nx = 4; Ny = 4; Nz = 2;
    }
    else if(ncube == 64) {
        Nx = 4; Ny = 4; Nz = 4;
    }

    hx = (Xlim[1]-Xlim[0])/Nx;
    hy = (Ylim[1]-Ylim[0])/Ny;
    hz = (Zlim[1]-Zlim[0])/Nz;

    std::vector<std::vector<AlfAtom>> splitAtoms;
    splitAtoms.resize(ncube);

    for (int j = 0; j < size; j++) {
        for (int k = 0; k < 3; k++) Point[k] = alfatoms[j].coord[k];
        Ix = (Point[0]-Xlim[0])/hx;
        Iy = (Point[1]-Ylim[0])/hy;
        Iz = (Point[2]-Zlim[0])/hz;
        idx = Ix + Iy*Nx + Iz*Nx*Ny;
        splitAtoms[idx].push_back(alfatoms[j]);
    }

    int nat = 0;
    Nval[0] = 0;
    for (int i = 0; i < ncube; i++) {
        Nval[i+1] = Nval[i] + splitAtoms[i].size();
        // std::cout << "i = " << i << " size = " << splitAtoms[i].size() << std::endl;
        for (int j = 0; j < splitAtoms[i].size(); j++) {
            alfatoms[nat+j] = splitAtoms[i][j];
        }
        nat += splitAtoms[i].size();
    }
}
}
