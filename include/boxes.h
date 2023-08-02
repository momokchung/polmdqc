///////////////////////////////////////////////////////////
//                                                       //
//  boxes.h  --  periodic boundary condition parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// xbox        length of a-axis of periodic box in Angstroms
// ybox        length of b-axis of periodic box in Angstroms
// zbox        length of c-axis of periodic box in Angstroms
// alphaA      angle between b- and c-axes of box in degrees
// betaA       angle between a- and c-axes of box in degrees
// gammaA      angle between a- and b-axes of box in degrees
// xbox2       half of the a-axis length of periodic box
// ybox2       half of the b-axis length of periodic box
// zbox2       half of the c-axis length of periodic box
// box34       three-fourths axis length of truncated octahedron
// volbox      volume in Ang**3 of the periodic box
// alpha_sin   sine of the alpha periodic box angle
// alpha_cos   cosine of the alpha periodic box angle
// beta_sin    sine of the beta periodic box angle
// beta_cos    cosine of the beta periodic box angle
// gamma_sin   sine of the gamma periodic box angle
// gamma_cos   cosine of the gamma periodic box angle
// beta_term   term used in generating triclinic box
// gamma_term  term used in generating triclinic box
// lvec        real space lattice vectors as matrix rows
// recip       reciprocal lattice vectors as matrix columns
// orthogonal  flag to mark periodic box as orthogonal
// monoclinic  flag to mark periodic box as monoclinic
// triclinic   flag to mark periodic box as triclinic
// octahedron  flag to mark box as truncated octahedron
// dodecadron  flag to mark box as rhombic dodecahedron
// nonprism    flag to mark octahedron or dodecahedron
// nosymm      flag to mark use or lack of lattice symmetry
// spacegrp    space group symbol for the unit cell type


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN double xbox,ybox,zbox;
QCMD_EXTERN double alphaA,betaA,gammaA;
QCMD_EXTERN double xbox2,ybox2,zbox2;
QCMD_EXTERN double box34,volbox;
QCMD_EXTERN double alpha_sin;
QCMD_EXTERN double alpha_cos;
QCMD_EXTERN double beta_sin;
QCMD_EXTERN double beta_cos;
QCMD_EXTERN double gamma_sin;
QCMD_EXTERN double gamma_cos;
QCMD_EXTERN double beta_term;
QCMD_EXTERN double gamma_term;
QCMD_EXTERN double lvec[3][3];
QCMD_EXTERN double recip[3][3];
QCMD_EXTERN bool orthogonal;
QCMD_EXTERN bool monoclinic;
QCMD_EXTERN bool triclinic;
QCMD_EXTERN bool octahedron;
QCMD_EXTERN bool dodecadron;
QCMD_EXTERN bool nonprism;
QCMD_EXTERN bool nosymm;
QCMD_EXTERN std::string spacegrp;
