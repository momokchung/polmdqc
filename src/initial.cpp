/////////////////////////////////////////////////////////
//                                                     //
//  initial.cpp  --  initial values and program setup  //
//                                                     //
/////////////////////////////////////////////////////////

// "initial" sets up original values for some parameters and
// variables that might not otherwise get initialized


#include "atoms.h"
#include "command.h"
#include "inform.h"
#include "initatom.h"
#include "initial.h"
#include "keys.h"
#include "params.h"
#include "promo.h"
#include "sizes.h"

void initial(int argc, char** argv)
{
    promo();

    command(argc, argv);

    initatom();

    // number of lines in the keyfile
    nkey = 0;

    // number of lines in the parameter file
    nprm = 0;

    // number of atoms in the system
    n = 0;

    // information levels within the program
    verbose = false;
    debug = false;
    silent = false;
    informAbort = false;
}
