// Author: Moses KJ Chung
// Year:   2024

#include "atomid.h"
#include "atoms.h"
#include "katoms.h"
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

    for (int i = 0; i < n; i++) {
        sym = symbol[i];
        atmn = symtoatmn[sym];
        atomic[i] = atmn;
        mass[i] = atmass[atmn-1];
    }
}
}
