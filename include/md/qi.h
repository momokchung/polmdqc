// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "realn.h"
#include <cmath>

namespace polmdqc
{
// Rt Q = G
inline void rotQI2GVector(const real (&rot)[3][3], real3 qif, real3& glf)
{
    glf = make_real3(dot3(rot[0][0], rot[1][0], rot[2][0], qif), dot3(rot[0][1], rot[1][1], rot[2][1], qif),
        dot3(rot[0][2], rot[1][2], rot[2][2], qif));
}

// R G = Q
void rotG2QIVector(const real (&rot)[3][3], real3 glf, real3& qif)
{
    qif = make_real3(dot3(rot[0][0], rot[0][1], rot[0][2], glf), dot3(rot[1][0], rot[1][1], rot[1][2], glf),
        dot3(rot[2][0], rot[2][1], rot[2][2], glf));
}

// R G Rt = Q
void rotG2QIMat_v1(const real (&rot)[3][3],
   real glxx, real glxy, real glxz, real glyy, real glyz, real glzz,
   real& qixx, real& qixy, real& qixz, real& qiyy, real& qiyz,
   real& qizz)
{
   real gl[3][3] = {{glxx, glxy, glxz}, {glxy, glyy, glyz}, {glxz, glyz, glzz}};
   real out[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
   // out[i][j] = sum(k,m) R[i][k] gl[k][m] Rt[m][j]
   //           = sum(k,m) R[i][k] gl[k][m] R[j][m]
   for (int i = 0; i < 2; ++i)
      for (int j = i; j < 3; ++j)
         for (int k = 0; k < 3; ++k)
            for (int m = 0; m < 3; ++m)
               out[i][j] += rot[i][k] * gl[k][m] * rot[j][m];
   qixx = out[0][0];
   qixy = out[0][1];
   qixz = out[0][2];
   qiyy = out[1][1];
   qiyz = out[1][2];
   // qizz = out[2][2];
   qizz = -(out[0][0] + out[1][1]);
}

// R G Rt = Q
void rotG2QIMat_v2(const real (&r)[3][3],
   real glxx, real glxy, real glxz, real glyy, real glyz, real glzz,
   real& qixx, real& qixy, real& qixz, real& qiyy, real& qiyz,
   real& qizz)
{
   // clang-format off
   qixx=r[0][0]*(r[0][0]*glxx+2*r[0][1]*glxy) + r[0][1]*(r[0][1]*glyy+2*r[0][2]*glyz) + r[0][2]*(r[0][2]*glzz+2*r[0][0]*glxz);
   qiyy=r[1][0]*(r[1][0]*glxx+2*r[1][1]*glxy) + r[1][1]*(r[1][1]*glyy+2*r[1][2]*glyz) + r[1][2]*(r[1][2]*glzz+2*r[1][0]*glxz);
   qixy=r[0][0]*(r[1][0]*glxx+r[1][1]*glxy+r[1][2]*glxz) + r[0][1]*(r[1][0]*glxy+r[1][1]*glyy+r[1][2]*glyz) + r[0][2]*(r[1][0]*glxz+r[1][1]*glyz+r[1][2]*glzz);
   qixz=r[0][0]*(r[2][0]*glxx+r[2][1]*glxy+r[2][2]*glxz) + r[0][1]*(r[2][0]*glxy+r[2][1]*glyy+r[2][2]*glyz) + r[0][2]*(r[2][0]*glxz+r[2][1]*glyz+r[2][2]*glzz);
   qiyz=r[1][0]*(r[2][0]*glxx+r[2][1]*glxy+r[2][2]*glxz) + r[1][1]*(r[2][0]*glxy+r[2][1]*glyy+r[2][2]*glyz) + r[1][2]*(r[2][0]*glxz+r[2][1]*glyz+r[2][2]*glzz);
   // clang-format on
   qizz = -(qixx + qiyy);
}

#define rotG2QIMatrix rotG2QIMat_v2
}