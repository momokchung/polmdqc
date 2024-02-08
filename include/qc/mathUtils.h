// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include <cmath>
#include <vector>

using std::sqrt;
using std::pow;
using std::exp;

namespace mathUtils{
/////////////////////////////////////////////////
//                                             //
//  mathUtils  --  Utility for math functions  //
//                                             //
/////////////////////////////////////////////////

void symmetrize(std::vector<std::vector<real>>& matrix);
void dsyevd(int N, real* A, real* w);
void dgemm(int m, int n, int k, real* A, real* B, real* C, char transa = 'N', char transb = 'N');
void triDiagSym(std::vector<std::vector<real>>& A);

////////////////////////////////////////////////////////////////
//                                                            //
//  inline int doubleFactorial  --  returns double factorial  //
//                                                            //
////////////////////////////////////////////////////////////////
inline int doubleFactorial(int n)
{
    int res = 1;
    for (int i = n; i >= 0; i = i - 2)
    {
        if (i == 0 || i == 1)
            return res;
        else
            res *= i;
    }
    return res;
}

//////////////////////////////////////////////////////////
//                                                      //
//  inline void multpoly  --  multiply two polynomials  //
//                                                      //
//////////////////////////////////////////////////////////

inline void multpoly(int iL, int jL, std::vector<real>& iPoly, std::vector<real>& jPoly, std::vector<real>& poly)
{
    for (int i = 0; i < iL; i++)
    {
        for (int j = 0; j < jL; j++)
        {
            poly[i + j] += iPoly[i] * jPoly[j];
        }
    }
    return;
}

///////////////////////////////////////////////////////
//                                                   //
//  inline void zero  --  make matrix elements zero  //
//                                                   //
///////////////////////////////////////////////////////

template <typename T>
inline void zero(int l, std::vector<T>& poly)
{
    for (int i = 0; i < l; i++)
    {
        poly[i] = 0;
    }
}

//////////////////////////////////////////////////////////
//                                                      //
//  inline real contract  --  contract two polynomials  //
//                                                      //
//////////////////////////////////////////////////////////

inline real contract(int l, std::vector<real>& poly1, std::vector<real>& poly2)
{
    real r = 0.;
    for (int i = 0; i < l; i++)
    {
        r += poly1[i] * poly2[i];
    }
    return r;
}

////////////////////////////////////////////////////////////
//                                                        //
//  inline void flatten  --  row major flatten 2D matrix  //
//                                                        //
////////////////////////////////////////////////////////////

inline void flatten(int m, int n, std::vector<std::vector<real>>& inputA, std::vector<real>& outputA)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            outputA[n * i + j] = inputA[i][j];
        }
    }
}
}
