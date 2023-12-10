// Author: Moses KJ Chung
// Year:   2023

#include "katomsqm.h"
#include "mathUtils.h"
#include "nuclearRepulsion.h"
#include <iostream>

namespace nuclearRepulsion
{
//////////////////////////////////////////////////////////////
//                                                          //
//  nuclearRepulsion  --  nuclear-nuclear repulsion energy  //
//                                                          //
//////////////////////////////////////////////////////////////

// nr    nuclear repulsion energy

real nr;

//////////////////////////////////////////////////////////////////////
//                                                                  //
//  nuclearRepulsion  --  compute nuclear-nuclear repulsion energy  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

void nuclearRepulsion()
{
    // get atom information
    int atomN = atoms::n;

    // outer loop over nucleus
    for (int i = 0; i < atomN; ++i)
    {
        int chargei = atoms::core[i];
        real coordxi = atoms::coordx[i];
        real coordyi = atoms::coordy[i];
        real coordzi = atoms::coordz[i];

        //inner loop over nucleus
        for (int j = 0; j < i; ++j)
        {
            int chargej = atoms::core[j];
            real coordxj = atoms::coordx[j];
            real coordyj = atoms::coordy[j];
            real coordzj = atoms::coordz[j];
            
            real xr = coordxj - coordxi;
            real yr = coordyj - coordyi;
            real zr = coordzj - coordzi;
            real r2 = xr * xr + yr * yr + zr * zr;
            real r = sqrt(r2);

            nr += chargei * chargej / r;
        }
    }

    // // print to debug
    // std::cout << "Nuclear Repulsion " << nr << std::endl;
}
}
