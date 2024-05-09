// Author: Moses KJ Chung
// Year:   2024

#include "inertia.h"
#include <algorithm>
#include <Eigen/Dense>

// remove
#include <iostream>

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  inertia  --  principal moments of inertia  //
//                                             //
/////////////////////////////////////////////////

// "inertia" computes the principal moments of inertia for the
// system, and optionally translates the center of mass to the
// origin and rotates the principal axes onto the global axes
//
//    mode = 1     print the moments and principal axes
//    mode = 2     move coordinates to standard orientation
//    mode = 3     perform both of the above operations
// 
// literature reference:
// 
// Herbert Goldstein, "Classical Mechanics, 2nd Edition",
// Addison-Wesley, Reading, MA, 1980; see the Euler angle
// xyz convention in Appendix B

void inertia(int mode, int n, real* mass, real* x, real* y, real* z)
{
    int index[3];
    real weigh,total,dot;
    real xcm,ycm,zcm;
    real xx,xy,xz,yy,yz,zz;
    real xterm,yterm,zterm;
    real phi,theta,psi;
    real moment[3],val[3],vec[3][3];
    Eigen::Matrix<real, 3, 3> evec,tensor;
    Eigen::Matrix<real, 3, 1> eval;

    // decide upon the type of output desired
    bool print = false;
    bool moved = false;
    if (mode==1 or mode==3) print = true;
    if (mode==2 or mode==3) moved = true;

    // compute the position of the center of mass
    total = 0;
    xcm = 0;
    ycm = 0;
    zcm = 0;
    for (int i = 0; i < n; i++) {
        weigh = mass[i];
        total += weigh;
        xcm += x[i]*weigh;
        ycm += y[i]*weigh;
        zcm += z[i]*weigh;
    }
    xcm /= total;
    ycm /= total;
    zcm /= total;

    // compute and then diagonalize the inertia tensor
    xx = 0;
    xy = 0;
    xz = 0;
    yy = 0;
    yz = 0;
    zz = 0;
    for (int i = 0; i < n; i++) {
        weigh = mass[i];
        xterm = x[i] - xcm;
        yterm = y[i] - ycm;
        zterm = z[i] - zcm;
        xx += xterm*xterm*weigh;
        xy += xterm*yterm*weigh;
        xz += xterm*zterm*weigh;
        yy += yterm*yterm*weigh;
        yz += yterm*zterm*weigh;
        zz += zterm*zterm*weigh;
    }
    tensor(0,0) = yy + zz;
    tensor(0,1) = -xy;
    tensor(0,2) = -xz;
    tensor(1,0) = -xy;
    tensor(1,1) = xx + zz;
    tensor(1,2) = -yz;
    tensor(2,0) = -xz;
    tensor(2,1) = -yz;
    tensor(2,2) = xx + yy;

    // eigenvalue problem
    Eigen::EigenSolver<Eigen::Matrix<real, 3, 3>> eigensolver;
    eigensolver.compute(tensor);
    eval = eigensolver.eigenvalues().real();
    evec = eigensolver.eigenvectors().real();

    // sort eigenvalue
    for (int i = 0; i < 3; i++) {
        index[i] = i;
        val[i] = eval(i,0);
    }
    std::sort(index, index + 3, [&](int i, int j) {
        return val[i] < val[j];
    });
    for (int i = 0; i < 3; i++) {
        int ind = index[i];
        moment[i] = val[ind];
        for (int j = 0; j < 3; j++) {
            vec[j][i] = evec(j,ind);
        }
    }

    // select the direction for each principal moment axis
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++) {
            xterm = vec[0][i] * (x[j]-xcm);
            yterm = vec[1][i] * (y[j]-ycm);
            zterm = vec[2][i] * (z[j]-zcm);
            dot = xterm + yterm + zterm;
            if (dot < 0.0) {
                for (int k = 0; k < 3; k++) {
                    vec[k][i] = -vec[k][i];
                }
            }
            if (dot != 0) break;
        }
    }

    // moment axes must give a right-handed coordinate system
    xterm = vec[0][0] * (vec[1][1]*vec[2][2]-vec[1][2]*vec[2][1]);
    yterm = vec[1][0] * (vec[0][2]*vec[2][1]-vec[0][1]*vec[2][2]);
    zterm = vec[2][0] * (vec[0][1]*vec[1][2]-vec[0][2]*vec[1][1]);
    dot = xterm + yterm + zterm;
    if (dot < 0) {
        for (int j = 0; j < 3; j++) {
            vec[j][2] = -vec[j][2];
        }
    }

    // principal moment axes form rows of Euler rotation matrix
    if (moved) {
        // translate to origin, then apply Euler rotation matrix
        for (int i = 0; i < n; i++) {
            xterm = x[i] - xcm;
            yterm = y[i] - ycm;
            zterm = z[i] - zcm;
            x[i] = vec[0][0]*xterm + vec[1][0]*yterm + vec[2][0]*zterm;
            y[i] = vec[0][1]*xterm + vec[1][1]*yterm + vec[2][1]*zterm;
            z[i] = vec[0][2]*xterm + vec[1][2]*yterm + vec[2][2]*zterm;
        }
    }
}
}
