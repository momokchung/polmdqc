//////////////////////////////////////////////////////////
//                                                      //
//  unitcell.cpp  --  get periodic boundary conditions  //
//                                                      //
//////////////////////////////////////////////////////////

// "unitcell" gets the periodic boundary box size and related
// values from an external keyword file


#include "bound.h"
#include "boxes.h"
#include "gettext.h"
#include "getword.h"
#include "keys.h"
#include "mathConst.h"
#include "unitcell.h"
#include "upcase.h"
#include <sstream>

void unitcell()
{
    // set the default values for periodic boundary conditions
    use_bounds = false;
    use_replica = false;

    // set the default values for the unit cell variables
    orthogonal = false;
    monoclinic = false;
    triclinic = false;
    octahedron = false;
    dodecadron = false;
    nonprism = false;
    nosymm = false;
    spacegrp = "          ";

    // get keywords containing crystal lattice dimensions
    for (int i = 0; i < nkey; i++) {
        std::string record = keyline[i];
        std::istringstream iss(record);
        std::string keyword;
        iss >> keyword;
        upcase(keyword);
        double nextDouble;
        if (keyword == "X-AXIS") {
            if (xbox == 0.) {
                if (iss >> nextDouble) xbox = nextDouble;
            }
        }
        else if (keyword == "Y-AXIS") {
            if (ybox == 0.)  {
                if (iss >> nextDouble) ybox = nextDouble;
            }
        }
        else if (keyword == "Z-AXIS") {
            if (zbox == 0.)  {
                if (iss >> nextDouble) zbox = nextDouble;
            }
        }
        else if (keyword == "A-AXIS") {
            if (xbox == 0.)  {
                if (iss >> nextDouble) xbox = nextDouble;
            }
        }
        else if (keyword == "B-AXIS") {
            if (ybox == 0.)  {
                if (iss >> nextDouble) ybox = nextDouble;
            }
        }
        else if (keyword == "C-AXIS") {
            if (zbox == 0.)  {
                if (iss >> nextDouble) zbox = nextDouble;
            }
        }
        else if (keyword == "ALPHA") {
            if (alpha == 0.)  {
                if (iss >> nextDouble) alpha = nextDouble;
            }
        }
        else if (keyword == "BETA") {
            if (beta == 0.)  {
                if (iss >> nextDouble) beta = nextDouble;
            }
        }
        else if (keyword == "GAMMA") {
            if (gamma == 0.)  {
                if (iss >> nextDouble) gamma = nextDouble;
            }
        }
        else if (keyword == "OCTAHEDRON") {
            octahedron = true;
        }
        else if (keyword == "DODECAHEDRON") {
            dodecadron = true;
        }
        else if (keyword == "NOSYMMETRY") {
            nosymm = true;
        }
        else if (keyword == "SPACEGROUP") {
            int next = 0;
            std::string nextString;
            iss >> nextString;
            getword(nextString, spacegrp, next);
        }
    }

    // use periodic boundary conditions if a cell was defined
    double boxmax = std::max({xbox,ybox,zbox});
    if (boxmax != 0.)  use_bounds = true;

    // set unspecified periodic boundary box lengths and angles
    if (use_bounds) {
        if (xbox !=  0.)  xbox = boxmax;
        if (ybox !=  0.)  ybox = boxmax;
        if (zbox !=  0.)  zbox = boxmax;
        if (alpha !=  0.)  alpha = 90.;
        if (beta !=  0.)  beta = 90.;
        if (gamma != 0.)  gamma = 90.;

        // determine the general periodic boundary lattice type
        if (nosymm) {
            triclinic = true;
        }
        else if (alpha==90. and beta==90. and gamma==90.) {
            orthogonal = true;
        }
        else if (alpha==90. and gamma==90.) {
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
        ybox = xbox;
        if (octahedron) zbox = xbox;
        else if (dodecadron) zbox = xbox * root2;
        alpha = 90.;
        beta = 90.;
        gamma = 90.;
    }
}
