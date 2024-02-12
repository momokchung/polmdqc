// Author: Moses KJ Chung
// Year:   2024

#include "atomid.h"
#include "atoms.h"
#include "ghost.h"
#include "ptable.h"

namespace polmdqc
{
//////////////////////////////////////////////
//                                          //
//  katomqm  --  atom parameter assignment  //
//                                          //
//////////////////////////////////////////////

// "katomqm" assigns atomic parameters to each atom

void katomqm()
{
    int atmn;
    std::string sym;

    // allocate global arrays from module atomid
    atomic.allocate(n);
    mass.allocate(n);

    // initialize atomic parameters
    for (int i = 0; i < n; i++) {
        atomic[i] = 0;
        mass[i] = 0.;
    }

    // assign atomic parameters
    for (int i = 0; i < n; i++) {
        sym = name[i];
        atmn = symtoatmn[sym];
        atomic[i] = atmn;
        if (!ghst[i]) {
            mass[i] = atmass[atmn-1];
        }
    }
}
}
