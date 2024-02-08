// Author: Moses KJ Chung
// Year:   2023

#include "align.h"
#include "atoms.h"
#include "bath.h"
#include "bound.h"
#include "boxes.h"
#include "cell.h"
#include "command.h"
#include "fft.h"
#include "files.h"
#include "group.h"
#include "inform.h"
#include "initatom.h"
#include "initres.h"
#include "keys.h"
#include "linmin.h"
#include "minima.h"
#include "molcul.h"
#include "mutant.h"
#include "neigh.h"
#include "openmp.h"
#include "output.h"
#include "params.h"
#include "pdb.h"
#include "promo.h"
#include "rigid.h"
#include "scales.h"
#include "sequen.h"
#include "socket.h"
#include "virial.h"
#include "warp.h"
#include "zclose.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  initial  --  initial values and program setup  //
//                                                 //
/////////////////////////////////////////////////////

// "initial" sets up original values for some parameters and
// variables that might not otherwise get initialized

void initial(int argc, char** argv)
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

    // names of biopolymer residue types
    initres();

    // number of lines in the keyfile
    nkey = 0;

    // number of lines in the parameter file
    nprm = 0;

    // number of atoms in the system
    n = 0;

    // number of molecules in the system
    nmol = 0;

    // number of unit cell and replicates
    ncell = 1;

    // number of atoms used in superposition
    nfit = 0;

    // number of mutated atoms in the system
    nmut = 0;

    // number of bonds added or deleted from Z-matrix
    nadd = 0;
    ndel = 0;

    // number of atoms in Protein Data Bank format
    npdb = 0;

    // number of residues and chains in biopolymer sequence
    nseq = 0;
    nchain = 0;

    // highest numbered previous cycle file
    nprior = 0;

    // pointer initialization for FFTW plans
    planf = 0;
    planb = 0;

    // information levels within the program
    verbose = false;
    debug = false;
    silent = false;
    informAbort = false;

    // integer flag for use of GPU coprocessor
    gpucard = 0;

    // flag for use of atom groups
    use_group = false;

    // flags for use of periodic boundaries
    use_bounds = false;
    use_replica = false;
    use_polymer = false;

    // flags for rebuilding of neighbor lists
    dovlst = true;
    dodlst = true;
    doclst = true;
    domlst = true;
    doulst = true;

    // flags for temperature and pressure baths
    isothermal = false;
    isobaric = false;

    // flag for use of internal virial
    use_virial = true;

    // flag for use of rigid bodies
    use_rigid = false;

    // flag to show setting of optimization scale factors
    set_scale = false;

    // flags for external Java socket communication
    sktstart = false;
    use_socket = false;

    // flags for potential energy smoothing
    use_smooth = false;
    use_dem = false;
    use_gda = false;
    use_tophat = false;
    use_stophat = false;

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
    stpmax = 0.;
}
}
