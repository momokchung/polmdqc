// Author: Moses KJ Chung
// Year:   2024

#include "atoms.h"
#include "basis.h"
#include "boxes.h"
#include "command.h"
#include "files.h"
#include "inform.h"
#include "initatom.h"
#include "initialqm.h"
#include "keys.h"
#include "minima.h"
#include "molcul.h"
#include "openmp.h"
#include "output.h"
#include "promo.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  initialqm  --  initial values and program setup  //
//                                                   //
///////////////////////////////////////////////////////

// "initialqm" sets up original values for some parameters
// and variables that might not otherwise get initialized

void initialqm(int argc, char** argv)
{
    // display program banner and copyright notice
    if (!test) promo();

    // command line arguments to the program
    command(argc, argv);

    // cores, thread count and options for OpenMP
    nproc = 1;
    nthread = 1;
    nproc = omp_get_num_procs();
    nthread = nproc;
    omp_set_num_threads(nthread);
    omp_set_nested(true);
#ifdef PolMDQC_ICPC
    // Intel compiler extensions to OpenMP standard, 268435456 bytes is 2**28 bytes, or 256 MB
    kmp_set_stacksize_s(268435456);
    kmp_set_blocktime(0);
#endif

    // atomic symbols, weights and radii
    initatom();

    // number of lines in the keyfile
    nkey = 0;

    // number of lines in the basis file
    nbss = 0;

    // number of atoms in the system
    n = 0;

    // number of molecules in the system
    nmol = 0;

    // highest numbered previous cycle file
    nprior = 0;

    // information levels within the program
    verbose = false;
    debug = false;
    silent = false;
    informAbort = false;

    // integer flag for use of GPU coprocessor
    gpucard = 0;

    // format for output of coordinates
    archive = true;
    binary = false;
    coordtype = "NONE";

    // default values for unit cell dimensions
    xbox = 0.;
    ybox = 0.;
    zbox = 0.;
    alphaA = 0.;
    betaA = 0.;
    gammaA = 0.;

    // default values used by optimizations
    fctmin = 0.;
    maxiter = 0;
    nextiter = 0;
    iprint = 0;
    iwrite = 0;
}
}
