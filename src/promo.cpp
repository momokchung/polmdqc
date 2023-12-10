// Author: Moses KJ Chung
// Year:   2023

#include "promo.h"
#include <stdio.h>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  promo  --  version info and copywrite notice  //
//                                                //
////////////////////////////////////////////////////

// "promo" writes a banner message containing information
// about the PolMDQC version, release date and copyright notice

void promo()
{
    printf("\n     ########################################################################");
    printf("\n   ###########################################################################");
    printf("\n  ###                                                                       ###");
    printf("\n ###               PolMDQC  ---  Software Tools for QC and MD                ###");
    printf("\n ##                                                                           ##");
    printf("\n ##                        Version 1.0.0    July 2023                         ##");
    printf("\n ##                                                                           ##");
    printf("\n ##          Copyright (c)  Moses KJ Chung & Jay W Ponder  2023-2023          ##");
    printf("\n ###                           All Rights Reserved                           ###");
    printf("\n  ###                                                                       ###");
    printf("\n   ###########################################################################");
    printf("\n     ########################################################################\n\n");
}
}
