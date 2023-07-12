///////////////////////////////////////////////////
//                                               //
//  torsions.cpp  --  locate and store torsions  //
//                                               //
///////////////////////////////////////////////////

// "torsions" finds the total number of torsional angles and
// the numbers of the four atoms defining each torsional angle


#include "atoms.h"
#include "bndstr.h"
#include "couple.h"
#include "fatal.h"
#include "tors.h"
#include "torsions.h"

void torsions()
{
    int maxtors;

    // perform dynamic allocation of some global arrays
    maxtors = 18 * n;
    itors.resize(maxtors,std::vector<int>(4));

    // loop over all bonds, storing the atoms in each torsion
    ntors = 0;
    for (int i = 0; i < nbond; i++) {
        int ib = ibnd[i][0];
        int ic = ibnd[i][1];
        for (int j = 0; j < n12[ib]; j++) {
            int ia = i12[ib][j];
            if (ia != ic) {
                for (int k = 0; k < n12[ic]; k++) {
                    int id = i12[ic][k];
                    if (id!=ib and id!=ia) {
                        ntors++;
                        if (ntors > maxtors) {
                            printf("\n TORSIONS  --  Too many Torsional Angles; Increase MAXTORS\n");
                            fatal();
                        }
                        itors[ntors-1][0] = ia;
                        itors[ntors-1][1] = ib;
                        itors[ntors-1][2] = ic;
                        itors[ntors-1][3] = id;
                    }
                }
            }
        }
    }
}
