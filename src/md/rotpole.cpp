// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "mathConst.h"
#include "mpole.h"
#include "repel.h"
#include "rotpole.h"
#include <cmath>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  rotpole  --  rotate multipoles to global frame  //
//                                                  //
//////////////////////////////////////////////////////

// "rotpole" constructs the global atomic multipoles by applying
// a rotation matrix to convert from local to global frame

void rotpole(RotMode rotMode)
{
    int i;
    Eigen::Matrix<real, 3, 3> a;
    bool planar;

    // rotate local multipoles to global frame at each site
    if (rotMode == RotMode::Mpole) {
        #pragma omp parallel for default(private)   \
        shared(n,pollist,pole,rpole)
        for (int i = 0; i < n; i++) {
            if (pollist[i] != -1) {
                rotmat(i,a,planar);
                rotsite(i,a,planar,pole,rpole);
            }
        }
    }
    else if (rotMode == RotMode::Repel) {
        #pragma omp parallel for default(private)   \
        shared(n,replist,repole,rrepole)
        for (int i = 0; i < n; i++) {
            if (replist[i] != -1) {
                rotmat(i,a,planar);
                rotsite(i,a,planar,repole,rrepole);
            }
        }
    }
}

//////////////////////////////////////////////////////
//                                                  //
//  rotrpole  --  rotate multipoles to local frame  //
//                                                  //
//////////////////////////////////////////////////////

// "rotrpole" constructs the local atomic multipoles by applying
// a rotation matrix to convert from global to local frame

void rotrpole(RotMode rotMode)
{
    int i;
    Eigen::Matrix<real, 3, 3> a;
    bool planar;

    // rotate global multipoles to local frame at each site
    if (rotMode == RotMode::Mpole) {
        #pragma omp parallel for default(private)   \
        shared(n,pollist,pole,rpole)
        for (int i = 0; i < n; i++) {
            if (pollist[i] != -1) {
                rotmat(i,a,planar);
                a = a.inverse().eval();
                planar = false;
                rotsite(i,a,planar,rpole,pole);
            }
        }
    }
    else if (rotMode == RotMode::Repel) {
        #pragma omp parallel for default(private)   \
        shared(n,replist,repole,rrepole)
        for (int i = 0; i < n; i++) {
            if (replist[i] != -1) {
                rotmat(i,a,planar);
                a = a.inverse().eval();
                planar = false;
                rotsite(i,a,planar,rrepole,repole);
            }
        }
    }
}

///////////////////////////////////////////////////
//                                               //
//  rotmat  --  local-to-global rotation matrix  //
//                                               //
///////////////////////////////////////////////////

// "rotmat" finds the rotation matrix that rotates the local
// coordinate system into the global frame at a specified atom

void rotmat(int i, Eigen::Matrix<real, 3, 3>& a, bool& planar)
{
    int ix,iy,iz;
    real r,dot;
    real eps,angle;
    real xi,yi,zi;
    real dx,dy,dz;
    real dx1,dy1,dz1;
    real dx2,dy2,dz2;
    real dx3,dy3,dz3;
    real dx4,dy4,dz4;
    LocalFrame axetyp;

    // get coordinates and frame definition for multipole site
    xi = x[i];
    yi = y[i];
    zi = z[i];
    iz = zaxis[i] - 1;
    ix = xaxis[i] - 1;
    iy = std::abs(yaxis[i]) - 1;
    axetyp = polaxe[i];
    planar = false;

    // use the identity matrix as the default rotation matrix
    a(0,0) = 1.;
    a(0,1) = 0.;
    a(0,2) = 0.;
    a(1,0) = 0.;
    a(1,1) = 1.;
    a(1,2) = 0.;
    a(2,0) = 0.;
    a(2,1) = 0.;
    a(2,2) = 1.;

    // get Z-Only rotation matrix elements for z-axis only
    if (axetyp == LocalFrame::ZOnly) {
        dx = x[iz] - xi;
        dy = y[iz] - yi;
        dz = z[iz] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(2,0) = dx / r;
        a(2,1) = dy / r;
        a(2,2) = dz / r;
        dx = 1.;
        dy = 0.;
        dz = 0.;
        dot = a(2,0);
        eps = 0.707;
        if (std::abs(dot) > eps) {
            dx = 0.;
            dy = 1.;
            dot = a(2,1);
        }
        dx = dx - dot*a(2,0);
        dy = dy - dot*a(2,1);
        dz = dz - dot*a(2,2);
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(0,0) = dx / r;
        a(0,1) = dy / r;
        a(0,2) = dz / r;
    }

    // get Z-then-X rotation matrix elements for z- and x-axes
    else if (axetyp == LocalFrame::ZthenX) {
        dx = x[iz] - xi;
        dy = y[iz] - yi;
        dz = z[iz] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(2,0) = dx / r;
        a(2,1) = dy / r;
        a(2,2) = dz / r;
        dx = x[ix] - xi;
        dy = y[ix] - yi;
        dz = z[ix] - zi;
        dot = dx*a(2,0) + dy*a(2,1) + dz*a(2,2);
        dx = dx - dot*a(2,0);
        dy = dy - dot*a(2,1);
        dz = dz - dot*a(2,2);
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(0,0) = dx / r;
        a(0,1) = dy / r;
        a(0,2) = dz / r;
    }

    // get Bisector rotation matrix elements for z- and x-axes
    else if (axetyp == LocalFrame::Bisector) {
        dx = x[iz] - xi;
        dy = y[iz] - yi;
        dz = z[iz] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx1 = dx / r;
        dy1 = dy / r;
        dz1 = dz / r;
        dx = x[ix] - xi;
        dy = y[ix] - yi;
        dz = z[ix] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx2 = dx / r;
        dy2 = dy / r;
        dz2 = dz / r;
        dx = dx1 + dx2;
        dy = dy1 + dy2;
        dz = dz1 + dz2;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(2,0) = dx / r;
        a(2,1) = dy / r;
        a(2,2) = dz / r;
        dot = dx2*a(2,0) + dy2*a(2,1) + dz2*a(2,2);
        dx = dx2 - dot*a(2,0);
        dy = dy2 - dot*a(2,1);
        dz = dz2 - dot*a(2,2);
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(0,0) = dx / r;
        a(0,1) = dy / r;
        a(0,2) = dz / r;
    }

    // get Z-Bisect rotation matrix elements for z- and x-axes;
    // use alternate x-axis if central atom is close to planar
    else if (axetyp == LocalFrame::ZBisect) {
        dx = x[iz] - xi;
        dy = y[iz] - yi;
        dz = z[iz] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(2,0) = dx / r;
        a(2,1) = dy / r;
        a(2,2) = dz / r;
        dx = x[ix] - xi;
        dy = y[ix] - yi;
        dz = z[ix] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx1 = dx / r;
        dy1 = dy / r;
        dz1 = dz / r;
        dx = x[iy] - xi;
        dy = y[iy] - yi;
        dz = z[iy] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx2 = dx / r;
        dy2 = dy / r;
        dz2 = dz / r;
        dx = dx1 + dx2;
        dy = dy1 + dy2;
        dz = dz1 + dz2;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx = dx / r;
        dy = dy / r;
        dz = dz / r;
        dot = dx*a(2,0) + dy*a(2,1) + dz*a(2,2);
        angle = 180. - radian*std::acos(dot);
        eps = 0.;
        if (angle < eps) {
            planar = true;
            dx = dy1*dz2 - dz1*dy2;
            dy = dz1*dx2 - dx1*dz2;
            dz = dx1*dy2 - dy1*dx2;
            dot = dx*a(2,0) + dy*a(2,1) + dz*a(2,2);
            if (dot < 0.) {
                dx = -dx;
                dy = -dy;
                dz = -dz;
                dot = -dot;
            }
        }
        dx = dx - dot*a(2,0);
        dy = dy - dot*a(2,1);
        dz = dz - dot*a(2,2);
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(0,0) = dx / r;
        a(0,1) = dy / r;
        a(0,2) = dz / r;
    }

    // get 3-Fold rotation matrix elements for z- and x-axes;
    // use alternate z-axis if central atom is close to planar
    else if (axetyp == LocalFrame::ThreeFold) {
        dx = x[iz] - xi;
        dy = y[iz] - yi;
        dz = z[iz] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx1 = dx / r;
        dy1 = dy / r;
        dz1 = dz / r;
        dx = x[ix] - xi;
        dy = y[ix] - yi;
        dz = z[ix] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx2 = dx / r;
        dy2 = dy / r;
        dz2 = dz / r;
        dx = x[iy] - xi;
        dy = y[iy] - yi;
        dz = z[iy] - zi;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        dx3 = dx / r;
        dy3 = dy / r;
        dz3 = dz / r;
        dx = dx1 + dx2 + dx3;
        dy = dy1 + dy2 + dy3;
        dz = dz1 + dz2 + dz3;
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        eps = 0.;
        if (r < eps) {
            planar = true;
            dx2 = x[ix] - x[iz];
            dy2 = y[ix] - y[iz];
            dz2 = z[ix] - z[iz];
            dx3 = x[iy] - x[iz];
            dy3 = y[iy] - y[iz];
            dz3 = z[iy] - z[iz];
            dx4 = dy2*dz3 - dz2*dy3;
            dy4 = dz2*dx3 - dx2*dz3;
            dz4 = dx2*dy3 - dy2*dx3;
            dot = dx4*dx + dy4*dy + dz4*dz;
            if (dot > 0.) {
                dx = dx4;
                dy = dy4;
                dz = dz4;
            }
            else {
                dx = -dx4;
                dy = -dy4;
                dz = -dz4;
            }
            r = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        a(2,0) = dx / r;
        a(2,1) = dy / r;
        a(2,2) = dz / r;
        dot = dx1*a(2,0) + dy1*a(2,1) + dz1*a(2,2);
        dx = dx1 - dot*a(2,0);
        dy = dy1 - dot*a(2,1);
        dz = dz1 - dot*a(2,2);
        r = std::sqrt(dx*dx + dy*dy + dz*dz);
        a(0,0) = dx / r;
        a(0,1) = dy / r;
        a(0,2) = dz / r;
    }

    // finally, find rotation matrix elements for the y-axis
    a(1,0) = a(0,2)*a(2,1) - a(0,1)*a(2,2);
    a(1,1) = a(0,0)*a(2,2) - a(0,2)*a(2,0);
    a(1,2) = a(0,1)*a(2,0) - a(0,0)*a(2,1);
}

/////////////////////////////////////////////////////
//                                                 //
//  rotsite  --  rotate input multipoles to final  //
//                                                 //
/////////////////////////////////////////////////////

// "rotsite" rotates atomic multipoles from the input to final
// frame at a specified atom by applying a rotation matrix

void rotsite(int ii, Eigen::Matrix<real, 3, 3>& a, bool& planar, MDQCArray2D<real,maxpole>& inpole, MDQCArray2D<real,maxpole>& outpole)
{
    int i,j,k,m;
    real spole[maxpole];
    real mp[3][3];
    real rp[3][3];
    LocalFrame axetyp;

    // copy input multipoles and modify at planar sites
    for (int i = 0; i < maxpole; i++) {
        spole[i] = inpole[ii][i];
    }
    if (planar) {
        axetyp = polaxe[ii];
        if (axetyp == LocalFrame::ZBisect) {
            spole[1] = 0.;
            spole[6] = 0.;
            spole[10] = 0.;
            spole[4] = 0.5 * (spole[4]+spole[8]);
            spole[8] = spole[4];
        }
        else if (axetyp == LocalFrame::ThreeFold) {
            for (int i = 1; i < maxpole; i++) {
                spole[i] = 0.;
            }
        }
    }

    // monopoles are the same in any coordinate frame
    outpole[ii][0] = spole[0];

    // rotate input dipoles to final coordinate frame
    for (int i = 1; i < 4; i++) {
        outpole[ii][i] = 0.;
        for (int j = 1; j < 4; j++) {
            outpole[ii][i] = outpole[ii][i] + spole[j]*a(j-1,i-1);
        }
    }

    // rotate input quadrupoles to final coordinate frame
    k = 4;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mp[j][i] = spole[k];
            rp[j][i] = 0.;
            k++;
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (j < i) {
                rp[j][i] = rp[i][j];
            }
            else {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        rp[j][i] = rp[j][i] + a(k,i)*a(m,j)*mp[m][k];
                    }
                }
            }
        }
    }
    k = 4;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            outpole[ii][k] = rp[j][i];
            k++;
        }
    }
}
}
