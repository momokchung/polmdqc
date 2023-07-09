/////////////////////////////////////////////////////////
//                                                     //
//  analyze.cpp  --  energy partitioning and analysis  //
//                                                     //
/////////////////////////////////////////////////////////

// "analyze" computes and displays the total potential energy;
// options are provided to display system and force field info,
// partition the energy by atom or by potential function type,
// show force field parameters by atom; output the large energy
// interactions and find electrostatic and inertial properties


#include "initial.h"
#include "getcart.h"
#include "mechanic.h"
#include <stdio.h>

int main(int argc, char** argv)
{
    initial(argc, argv);

    int ixyz;
    getcart(ixyz);
    mechanic();

    printf("\n");
}
