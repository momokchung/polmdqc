////////////////////////////////////////////////////////
//                                                    //
//  mechanic.cpp  --  initialize molecular mechanics  //
//                                                    //
////////////////////////////////////////////////////////

// "mechanic" sets up needed parameters for the potential energy
// calculation and reads in many of the user selectable options


#include "active.h"
#include "angles.h"
#include "attach.h"
#include "bitors.h"
#include "bonds.h"
#include "fatal.h"
#include "inform.h"
#include "lattice.h"
#include "mechanic.h"
#include "torsions.h"
#include "unitcell.h"

void mechanic()
{
    // set the bonded connectivity lists and active atoms
    attach();
    active();

    // find bonds, angles, torsions, bitorsions and small rings
    bonds();
    angles();
    torsions();
    bitors();
    // rings();

    // get the base force field from parameter file and keyfile
    // field(); do

    // find unit cell type, lattice parameters and cutoff values
    unitcell();
    lattice();
    // polymer();
    // cutoffs(); do

    // setup needed for potential energy smoothing methods
    // flatten();

    // assign atom types, classes and other atomic information
    // katom(); do

    // assign atoms to molecules and set the atom groups
    // molecule(); do
    // cluster(); do

    // find any pisystem atoms, bonds and torsional angles
    // orbital();

    // assign bond, angle and cross term potential parameters
    // kbond();
    // kangle();
    // kstrbnd();
    // kurey();
    // kangang();

    // assign out-of-plane deformation potential parameters
    // kopbend();
    // kopdist();
    // kimprop();
    // kimptor();

    // assign torsion and torsion cross term potential parameters
    // ktors();
    // kpitors();
    // kstrtor();
    // kangtor();
    // ktortor();

    // assign electrostatic interaction potential parameters
    // kcharge();
    // kdipole();
    // kmpole(); do
    // kpolar(); do
    // kchgtrn();
    // kchgflx();

    // assign van der Waals, repulsion and dispersion parameters
    // kvdw(); do
    // krepel();
    // kdisp();

    // assign solvation, metal, pisystem and restraint parameters
    // ksolv(); do
    // kmetal();
    // korbit();
    // kgeom();
    // kextra();

    // assign electrostatic and dispersion Ewald sum parameters
    // kewald(); do

    // set any holonomic interatomic distance constraints
    // shakeup();

    // set hybrid parameter values for free energy perturbation
    // mutate();

    // quit if essential parameter information is missing
    if (informAbort) {
        printf("\n MECHANIC  --  Some Required Potential Energy Parameters are Undefined\n");
        fatal();
    }
}
