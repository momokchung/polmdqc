// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "damping.h"
#include "libfunc.h"
#include "mplpot.h"
#include "qi.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////////
//                                                          //
//  pairMpoleQI  --  pairwise QI electrostatic calculation  //
//                                                          //
//////////////////////////////////////////////////////////////

// "pairMpoleQI" calculates the pairwise electrostatic energy 
// and/or gradient due to atomic multipole interactions using
// quasi-internal frames

template <bool do_e, bool do_g, bool do_v>
inline void pairMpoleQI(
    real r2, real xr, real yr, real zr, real mscale,
    real ci, real Idx, real Idy, real Idz, real Iqxx, real Iqxy, real Iqxz,
    real Iqyy, real Iqyz, real Iqzz,
    real ck, real Kdx, real Kdy, real Kdz, real Kqxx, real Kqxy, real Kqxz,
    real Kqyy, real Kqyz, real Kqzz,
    real f,
    real& frcx, real& frcy, real& frcz,
    real& trqxi, real& trqyi, real& trqzi,
    real& trqxk, real& trqyk, real& trqzk,
    real& e,
    real& vxx, real& vxy, real& vxz,
    real& vyy, real& vyz, real& vzz)
{
    real3 dR = make_real3(xr, yr, zr);
    real3 Id = make_real3(Idx, Idy, Idz);
    real3 Kd = make_real3(Kdx, Kdy, Kdz);

    // a rotation matrix that rotates (xr,yr,zr) to (0,0,r); R G = Q
    real rot[3][3];
    real bn[6];
    real sr3, sr5, sr7, sr9;
    real r = REAL_SQRT(r2);
    real invr1 = 1/r;
    real rr2 = invr1 * invr1;
    real rr1 = mscale * f * invr1;
    real rr3 = rr1 * rr2;
    real rr5 = 3 * rr3 * rr2;
    real rr7 = 5 * rr5 * rr2;
    real rr9 = 7 * rr7 * rr2;
    real rr11;
    if constexpr (do_g) rr11 = 9 * rr9 * rr2;

    real3 rotz = invr1 * dR;
    // pick a random vector as rotx; rotx and rotz cannot be parallel
    real3 rotx = rotz;
    if (dR.y != 0 || dR.z != 0)
        rotx.x += 1;
    else
        rotx.y += 1;
    // Gram–Schmidt process for rotx with respect to rotz
    rotx -= dot3(rotx, rotz) * rotz;
    // normalize rotx
    real invxlen = 1/REAL_SQRT(dot3(rotx, rotx));
    rotx = invxlen * rotx;
    real3 roty = cross(rotz, rotx);
    rot[0][0] = rotx.x;
    rot[0][1] = rotx.y;
    rot[0][2] = rotx.z;
    rot[1][0] = roty.x;
    rot[1][1] = roty.y;
    rot[1][2] = roty.z;
    rot[2][0] = rotz.x;
    rot[2][1] = rotz.y;
    rot[2][2] = rotz.z;

    real3 di, dk;
    rotG2QIVector(rot, Id, di);
    rotG2QIVector(rot, Kd, dk);
    real qixx, qixy, qixz, qiyy, qiyz, qizz;
    real qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;
    rotG2QIMatrix(rot, Iqxx, Iqxy, Iqxz, Iqyy, Iqyz, Iqzz, qixx, qixy, qixz, qiyy, qiyz, qizz);
    rotG2QIMatrix(rot, Kqxx, Kqxy, Kqxz, Kqyy, Kqyz, Kqzz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz);

    // phi,dphi/d(x,y,z),d2phi/dd(xx,yy,zz,xy,xz,yz)
    //   0        1 2 3            4  5  6  7  8  9
    real phi1[10] = {0};
    real phi2[10] = {0};
    real phi1z[10] = {0};

    // C-C
    {
        real coef1 = rr1;
        real coef3 = rr3 * r;
        // phi_c c
        phi1[0] += coef1 * ck;
        phi2[0] += coef1 * ci;
        phi1z[0] += coef3 * ck;
    }

    // D-C and C-D
    {
        real coef3 = rr3 * r;
        real coef5 = (rr3 - rr5 * r2);
        // phi_d c
        phi1[0] += -coef3 * dk.z;
        phi2[0] += coef3 * di.z;
        phi1z[0] += coef5 * dk.z;
        // dphi_c d
        // phi1[1]; phi1[2];
        phi1[3] += coef3 * ck;
        // phi2[1]; phi2[2];
        phi2[3] += -coef3 * ci;
        // phi1z[1]; phi1z[2];
        phi1z[3] += -coef5 * ck;
    }

    // D-D
    {
        real coef3 = rr3;
        real coef5 = (rr3 - rr5 * r2);
        real coez5 = rr5 * r;
        real coez7 = (3 * rr5 - rr7 * r2) * r;
        // dphi_d d
        phi1[1] += coef3 * dk.x;
        phi1[2] += coef3 * dk.y;
        phi1[3] += coef5 * dk.z;
        phi2[1] += coef3 * di.x;
        phi2[2] += coef3 * di.y;
        phi2[3] += coef5 * di.z;
        phi1z[1] += coez5 * dk.x;
        phi1z[2] += coez5 * dk.y;
        phi1z[3] += coez7 * dk.z;
    }

    // Q-C and C-Q
    {
        real coef5 = rr5 * r2;
        real coez5 = rr5 * r;
        real coez7 = rr7 * r2 * r;
        // phi_q c
        phi1[0] += coef5 * qkzz;
        phi2[0] += coef5 * qizz;
        phi1z[0] += -(2 * coez5 - coez7) * qkzz;
        // d2phi_c q
        // phi1[4]; phi1[5];
        phi1[6] += coef5 * ck;
        // phi1[7]; phi1[8]; phi1[9];
        // phi2[4]; phi2[5];
        phi2[6] += coef5 * ci;
        // phi2[7]; phi2[8]; phi2[9];
        // phi1z[4]; phi1z[5];
        phi1z[6] += -(2 * coez5 - coez7) * ck;
        // phi1z[7]; phi1z[8]; phi1z[9];
    }

    // Q-D and D-Q
    {
        real coef5 = rr5 * r;
        real coef7 = rr7 * r2 * r;
        real coez7 = (rr5 - rr7 * r2);
        real coez9 = (3 * rr7 - rr9 * r2) * r2;
        // dphi_q d
        phi1[1] += -2 * coef5 * qkxz;
        phi1[2] += -2 * coef5 * qkyz;
        phi1[3] += -(2 * coef5 - coef7) * qkzz;
        phi2[1] += 2 * coef5 * qixz;
        phi2[2] += 2 * coef5 * qiyz;
        phi2[3] += (2 * coef5 - coef7) * qizz;
        phi1z[1] += 2 * coez7 * qkxz;
        phi1z[2] += 2 * coez7 * qkyz;
        phi1z[3] += (2 * coez7 - coez9) * qkzz;
        // d2phi_d q
        // phi1[4]; phi1[5];
        phi1[6] += (2 * coef5 - coef7) * dk.z;
        // phi1[7];
        phi1[8] += 2 * coef5 * dk.x;
        phi1[9] += 2 * coef5 * dk.y;
        //
        // phi2[4]; phi2[5];
        phi2[6] += -(2 * coef5 - coef7) * di.z;
        // phi2[7];
        phi2[8] += -2 * coef5 * di.x;
        phi2[9] += -2 * coef5 * di.y;
        //
        // phi1z[4]; phi1z[5];
        phi1z[6] += -(2 * coez7 - coez9) * dk.z;
        // phi1z[7];
        phi1z[8] += -2 * coez7 * dk.x;
        phi1z[9] += -2 * coez7 * dk.y;
    }

    // Q-Q
    {
        // d2phi_q q
        real coef5 = rr5;
        real coef7 = rr7 * r2;
        real coef9 = rr9 * r2 * r2;
        real coez7 = rr7 * r;
        real coez9 = rr9 * r2 * r;
        real coez11 = rr11 * r2 * r2 * r;
        //
        phi1[4] += 2 * coef5 * qkxx;
        phi1[5] += 2 * coef5 * qkyy;
        phi1[6] += (2 * coef5 - 4 * coef7 + coef9) * qkzz;
        phi1[7] += 4 * coef5 * qkxy;
        phi1[8] += 4 * (coef5 - coef7) * qkxz;
        phi1[9] += 4 * (coef5 - coef7) * qkyz;
        //
        phi2[4] += 2 * coef5 * qixx;
        phi2[5] += 2 * coef5 * qiyy;
        phi2[6] += (2 * coef5 - 4 * coef7 + coef9) * qizz;
        phi2[7] += 4 * coef5 * qixy;
        phi2[8] += 4 * (coef5 - coef7) * qixz;
        phi2[9] += 4 * (coef5 - coef7) * qiyz;
        //
        phi1z[4] += 2 * coez7 * qkxx;
        phi1z[5] += 2 * coez7 * qkyy;
        phi1z[6] += (10 * coez7 - 8 * coez9 + coez11) * qkzz;
        phi1z[7] += 4 * coez7 * qkxy;
        phi1z[8] += 4 * (3 * coez7 - coez9) * qkxz;
        phi1z[9] += 4 * (3 * coez7 - coez9) * qkyz;
    }

    if constexpr (do_e) {
        e = phi1[0] * ci + phi1[1] * di.x + phi1[2] * di.y + phi1[3] * di.z + phi1[4] * qixx + phi1[5] * qiyy
            + phi1[6] * qizz + phi1[7] * qixy + phi1[8] * qixz + phi1[9] * qiyz;
    }

    real3 frc, trq1, trq2;
    if constexpr (do_g) {
        // torque
        real3 trqa = cross(phi1[1], phi1[2], phi1[3], di);
        trqa.x += phi1[9] * (qizz - qiyy) + 2 * (phi1[5] - phi1[6]) * qiyz + phi1[7] * qixz - phi1[8] * qixy;
        trqa.y += phi1[8] * (qixx - qizz) + 2 * (phi1[6] - phi1[4]) * qixz + phi1[9] * qixy - phi1[7] * qiyz;
        trqa.z += phi1[7] * (qiyy - qixx) + 2 * (phi1[4] - phi1[5]) * qixy + phi1[8] * qiyz - phi1[9] * qixz;
        real3 trqb = cross(phi2[1], phi2[2], phi2[3], dk);
        trqb.x += phi2[9] * (qkzz - qkyy) + 2 * (phi2[5] - phi2[6]) * qkyz + phi2[7] * qkxz - phi2[8] * qkxy;
        trqb.y += phi2[8] * (qkxx - qkzz) + 2 * (phi2[6] - phi2[4]) * qkxz + phi2[9] * qkxy - phi2[7] * qkyz;
        trqb.z += phi2[7] * (qkyy - qkxx) + 2 * (phi2[4] - phi2[5]) * qkxy + phi2[8] * qkyz - phi2[9] * qkxz;
        trq1 = trqa;
        trq2 = trqb;

        // gradient
        real frc1z = phi1z[0] * ci + phi1z[1] * di.x + phi1z[2] * di.y + phi1z[3] * di.z + phi1z[4] * qixx
            + phi1z[5] * qiyy + phi1z[6] * qizz + phi1z[7] * qixy + phi1z[8] * qixz + phi1z[9] * qiyz;
        frc.x = -invr1 * (trqa.y + trqb.y);
        frc.y = invr1 * (trqa.x + trqb.x);
        frc.z = frc1z;

        real3 glfrc;
        rotQI2GVector(rot, frc, glfrc);
        frc = glfrc;
        frcx = frc.x;
        frcy = frc.y;
        frcz = frc.z;
        real3 gltrq1;
        rotQI2GVector(rot, trq1, gltrq1);
        trqxi = gltrq1.x;
        trqyi = gltrq1.y;
        trqzi = gltrq1.z;
        real3 gltrq2;
        rotQI2GVector(rot, trq2, gltrq2);
        trqxk = gltrq2.x;
        trqyk = gltrq2.y;
        trqzk = gltrq2.z;
    }
    if constexpr (do_v) {
        vxx = -dR.x * frc.x;
        vxy = (real)-0.5 * (dR.y * frc.x + dR.x * frc.y);
        vxz = (real)-0.5 * (dR.z * frc.x + dR.x * frc.z);
        vyy = -dR.y * frc.y;
        vyz = (real)-0.5 * (dR.z * frc.y + dR.y * frc.z);
        vzz = -dR.z * frc.z;
    }
}
}
