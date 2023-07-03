/////////////////////////////////////////////////////////
//                                                     //
//  nuclear.h  --  nuclear electron attraction matrix  //
//                                                     //
/////////////////////////////////////////////////////////


#pragma once
#include "init.h"
#include <vector>

namespace nuclear
{
extern std::vector< std::vector<real> > cartNE;
extern std::vector< std::vector<real> > sphNE;

void nuclearOS();


///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//  void recursionNeA  --  first recursion from Obara Saika nuclear electron attraction  //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

inline void recursionNeA(int i, int j, real xPA, real xPC, real aP2, std::vector< std::vector< std::vector<real> > >& xNe)
{
    if (i == 0) return;
    
    int maxt = i + j;

    for (int t = 0; t < maxt; ++t)
    {
        real tmp = xNe[i - 1][j][t];
        xNe[i][j][t] += xPA * tmp;
        xNe[i][j][t + 1] -= xPC * tmp;
    }
    if (i > 1)
    {
        for (int t = 0; t < maxt; ++t)
        {
            real tmp = (i - 1) / aP2 * xNe[i - 2][j][t];
            xNe[i][j][t] += tmp;
            xNe[i][j][t + 1] -= tmp;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//  void recursionNeB  --  second recursion from Obara Saika nuclear electron attraction  //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

inline void recursionNeB(int i, int j, real xPB, real xPC, real aP2, std::vector< std::vector< std::vector<real> > >& xNe)
{
    if (j == 0) return;

    int maxt = i + j;

    for (int t = 0; t < maxt; ++t)
    {
        real tmp = xNe[i][j - 1][t];
        xNe[i][j][t] += xPB * tmp;
        xNe[i][j][t + 1] -= xPC * tmp;
    }
    if (j > 1)
    {
        for (int t = 0; t < maxt; ++t)
        {
            real tmp = (j - 1) / aP2 * xNe[i][j - 2][t];
            xNe[i][j][t] += tmp;
            xNe[i][j][t + 1] -= tmp;
        }
    }
    if (i > 0)
    {
        for (int t = 0; t < maxt; ++t)
        {
            real tmp = i / aP2 * xNe[i - 1][j - 1][t];
            xNe[i][j][t] += tmp;
            xNe[i][j][t + 1] -= tmp;
        }
    }
}
}
