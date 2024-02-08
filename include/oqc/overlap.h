// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include <vector>

namespace overlap
{
///////////////////////////////////
//                               //
//  overlap  --  overlap matrix  //
//                               //
///////////////////////////////////

extern std::vector<std::vector<real>> cartS;
extern std::vector<std::vector<real>> sphS;

void overlapOS();

//////////////////////////////////////////////////////////////////////
//                                                                  //
//  void recursionOA  --  first recursion from Obara Saika overlap  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

inline void recursionOA(int i, int j, real xPA, real aP2, std::vector<std::vector<real>>& sx)
{
    if (i == 0) return;

    sx[i][j] += xPA * sx[i - 1][j];
    if (i > 1)
    {
        sx[i][j] += (i - 1) / aP2 * sx[i - 2][j];
    }
    return;
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  void recursionOB  --  second recursion from Obara Saika overlap  //
//                                                                   //
///////////////////////////////////////////////////////////////////////

inline void recursionOB(int i, int j, real xPB, real aP2, std::vector<std::vector<real>>& sx)
{
    if (j == 0) return;

    sx[i][j] += xPB * sx[i][j - 1];
    if (j > 1)
    {
        sx[i][j] += (j - 1) / aP2 * sx[i][j - 2];
    }
    if (i > 0)
    {
        sx[i][j] += i / aP2 * sx[i - 1][j - 1];
    }
    return;
}
}
