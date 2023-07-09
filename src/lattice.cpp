///////////////////////////////////////////////////////////
//                                                       //
//  lattice.cpp  --  setup periodic boundary conditions  //
//                                                       //
///////////////////////////////////////////////////////////

// "lattice" stores the periodic box dimensions and sets angle
// values to be used in computing fractional coordinates


#include "bound.h"
#include "boxes.h"
#include "cell.h"
#include "inform.h"
#include "lattice.h"
#include "mathConst.h"
#include <cmath>

void lattice()
{
    // use periodic boundary conditions if a cell was defined
    double boxmax = std::max({xbox, ybox, zbox});
    if (boxmax != 0.0) use_bounds = true;

    // set unspecified periodic boundary box lengths and angles
    if (use_bounds) {
        if (xbox == 0.)  xbox = boxmax;
        if (ybox == 0.)  ybox = boxmax;
        if (zbox == 0.)  zbox = boxmax;
        if (alpha == 0.)  alpha = 90.;
        if (beta == 0.)  beta = 90.;
        if (gamma == 0.)  gamma = 90.;

        // determine the general periodic boundary lattice type
        orthogonal = false;
        monoclinic = false;
        triclinic = false;
        if (nosymm) {
            triclinic = true;
        }
        else if (alpha == 90. && beta == 90. && gamma == 90.) {
            orthogonal = true;
        }
        else if (alpha == 90.0 && gamma == 90.0) {
            monoclinic = true;
        }
        else {
            triclinic = true;
        }
    }

    // set lattice values for non-prism periodic boundaries
    if (octahedron or dodecadron) {
        orthogonal = false;
        monoclinic = false;
        triclinic = false;
        nonprism = true;
    }

    // compute and store half box lengths and other lengths
    xbox2 = 0.5 * xbox;
    ybox2 = 0.5 * ybox;
    zbox2 = 0.5 * zbox;
    if (octahedron)  box34 = 0.75 * xbox;

    // set replicated cell dimensions equal to the unit cell
    xcell = xbox;
    ycell = ybox;
    zcell = zbox;
    xcell2 = xbox2;
    ycell2 = ybox2;
    zcell2 = zbox2;

    // get values needed for fractional coordinate computations
    if (triclinic) {
        alpha_sin = std::sin(alpha/radian);
        alpha_cos = std::cos(alpha/radian);
        beta_sin = std::sin(beta/radian);
        beta_cos = std::cos(beta/radian);
        gamma_sin = std::sin(gamma/radian);
        gamma_cos = std::cos(gamma/radian);
        beta_term = (alpha_cos - beta_cos*gamma_cos) / gamma_sin;
        gamma_term = std::sqrt(std::pow(beta_sin,2) - std::pow(beta_term,2));
    }
    else if (monoclinic) {
        alpha_sin = 1.;
        alpha_cos = 0.;
        beta_sin = std::sin(beta/radian);
        beta_cos = std::cos(beta/radian);
        gamma_sin = 1.;
        gamma_cos = 0.;
        beta_term = 0.;
        gamma_term = beta_sin;
    }
    else {
        alpha_sin = 1.;
        alpha_cos = 0.;
        beta_sin = 1.;
        beta_cos = 0.;
        gamma_sin = 1.;
        gamma_cos = 0.;
        beta_term = 0.;
        gamma_term = 1.;
    }

    // determine the volume of the parent periodic box
    volbox = 0.;
    if (triclinic) {
        volbox = (gamma_sin*gamma_term) * xbox * ybox * zbox;
    }
    else if (monoclinic) {
        volbox = beta_sin * xbox * ybox * zbox;
    }
    else {
        volbox = xbox * ybox * zbox;
    }

    // compute and store real space lattice vectors as rows
    double ar1 = xbox;
    double ar2 = 0.;
    double ar3 = 0.;
    double br1 = ybox * gamma_cos;
    double br2 = ybox * gamma_sin;
    double br3 = 0.;
    double cr1 = zbox * beta_cos;
    double cr2 = zbox * beta_term;
    double cr3 = zbox * gamma_term;
    lvec[0][0] = ar1;
    lvec[1][0] = ar2;
    lvec[2][0] = ar3;
    lvec[0][1] = br1;
    lvec[1][1] = br2;
    lvec[2][1] = br3;
    lvec[0][2] = cr1;
    lvec[1][2] = cr2;
    lvec[2][2] = cr3;

    // compute and store reciprocal lattice vectors as columns
    if (volbox != 0.) {
        recip[0][0] = (br2*cr3 - cr2*br3) / volbox;
        recip[0][1] = (br3*cr1 - cr3*br1) / volbox;
        recip[0][2] = (br1*cr2 - cr1*br2) / volbox;
        recip[1][0] = (cr2*ar3 - ar2*cr3) / volbox;
        recip[1][1] = (cr3*ar1 - ar3*cr1) / volbox;
        recip[1][2] = (cr1*ar2 - ar1*cr2) / volbox;
        recip[2][0] = (ar2*br3 - br2*ar3) / volbox;
        recip[2][1] = (ar3*br1 - br3*ar1) / volbox;
        recip[2][2] = (ar1*br2 - br1*ar2) / volbox;
    }

    // correct volume of non-parallelepiped periodic cells
    if (nonprism)  volbox = 0.5 * volbox;
}
