// Author: Moses KJ Chung
// Year:   2023

#pragma once

namespace polmdqc
{
//////////////////////////////////////////////
//                                          //
//  damping  --  various damping functions  //
//                                          //
//////////////////////////////////////////////

void damppole(double r,int rorder, double alphai, double alphak, double* dmpi, double* dmpk, double* dmpik);
}
