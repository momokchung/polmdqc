// Author: Moses KJ Chung
// Year:   2023

#include "bound.h"
#include "boxes.h"
#include "gettext.h"
#include "getword.h"
#include "keys.h"
#include "mathConst.h"
#include "unitcell.h"
#include "upcase.h"
#include <algorithm>
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  unitcell  --  get periodic boundary conditions  //
//                                                  //
//////////////////////////////////////////////////////

// "unitcell" gets the periodic boundary box size and related
// values from an external keyword file

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
        std::string keyword = "";
        iss >> keyword;
        upcase(keyword);
        real nextDouble;
        if (keyword == "X-AXIS") {
            if (xbox == 0.) {
                if (iss >> nextDouble) xbox = nextDouble;
            }
        }
        else if (keyword == "Y-AXIS") {
            if (ybox == 0.) {
                if (iss >> nextDouble) ybox = nextDouble;
            }
        }
        else if (keyword == "Z-AXIS") {
            if (zbox == 0.) {
                if (iss >> nextDouble) zbox = nextDouble;
            }
        }
        else if (keyword == "A-AXIS") {
            if (xbox == 0.) {
                if (iss >> nextDouble) xbox = nextDouble;
            }
        }
        else if (keyword == "B-AXIS") {
            if (ybox == 0.) {
                if (iss >> nextDouble) ybox = nextDouble;
            }
        }
        else if (keyword == "C-AXIS") {
            if (zbox == 0.) {
                if (iss >> nextDouble) zbox = nextDouble;
            }
        }
        else if (keyword == "ALPHA") {
            if (alphaA == 0.) {
                if (iss >> nextDouble) alphaA = nextDouble;
            }
        }
        else if (keyword == "BETA") {
            if (betaA == 0.) {
                if (iss >> nextDouble) betaA = nextDouble;
            }
        }
        else if (keyword == "GAMMA") {
            if (gammaA == 0.) {
                if (iss >> nextDouble) gammaA = nextDouble;
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
            std::string nextString = "";
            iss >> nextString;
            getword(nextString, spacegrp, next);
        }
    }

    // use periodic boundary conditions if a cell was defined
    real boxmax = std::max({xbox, ybox, zbox});
    if (boxmax != 0.) use_bounds = true;

    // set unspecified periodic boundary box lengths and angles
    if (use_bounds) {
        if (xbox == 0.) xbox = boxmax;
        if (ybox == 0.) ybox = boxmax;
        if (zbox == 0.) zbox = boxmax;
        if (alphaA == 0.) alphaA = 90.;
        if (betaA == 0.) betaA = 90.;
        if (gammaA == 0.) gammaA = 90.;

        // determine the general periodic boundary lattice type
        if (nosymm) {
            triclinic = true;
        }
        else if (alphaA == 90. and betaA == 90. and gammaA == 90.) {
            orthogonal = true;
        }
        else if (alphaA == 90. and gammaA == 90.) {
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
        alphaA = 90.;
        betaA = 90.;
        gammaA = 90.;
    }
}
}
