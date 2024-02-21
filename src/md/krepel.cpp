// Author: Moses KJ Chung
// Year:   2024

#include "atomid.h"
#include "atoms.h"
#include "chkpole.h"
#include "gettext.h"
#include "inform.h"
#include "krepel.h"
#include "krepl.h"
#include "keys.h"
#include "mpole.h"
#include "potent.h"
#include "repel.h"
#include "reppot.h"
#include "sizes.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  krepel  --  Pauli repulsion term assignment  //
//                                               //
///////////////////////////////////////////////////

// "krepel" assigns the size values, exponential parameter and
// number of valence electrons for Pauli repulsion interactions
// and processes any new or changed values for these parameters

void krepel()
{
    int k;
    int ia,ic,next;
    real spr,apr,epr;
    bool header;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // process keywords containing Pauli repulsion parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "REPULSION") {
            k = 0;
            spr = 0.;
            apr = 0.;
            epr = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> k >> spr >> apr >> epr;
            if (k > 0) {
                k--;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Pauli Repulsion Parameters :\n\n");
                    printf("     Atom Class               Size           Damp        Valence\n\n");
                }
                if (k < maxclass) {
                    prsiz[k] = spr;
                    prdmp[k] = apr;
                    prele[k] = -std::abs(epr);
                    if (!silent) {
                        printf("      %6d       %15.4f%15.4f%15.3f\n", k+1, spr, apr, epr);
                    }
                }
                else {
                    printf("\n KREPEL  --  Too many Pauli Repulsion Parameters\n");
                    informAbort = true;
                }
            }
        }
    }

    // perform dynamic allocation of some global arrays
    irep.allocate(n);
    replist.allocate(n);
    sizpr.allocate(n);
    dmppr.allocate(n);
    elepr.allocate(n);
    repole.allocate(n);
    rrepole.allocate(n);

    // assign the repulsion size, alpha and valence parameters
    for (int i = 0; i < n; i++) {
        irep[i] = -1;
        replist[i] = -1;
        sizpr[i] = 0.;
        dmppr[i] = 0.;
        elepr[i] = 0.;
        ic = atomClass[i];
        if (ic != -1) {
            sizpr[i] = prsiz[ic];
            dmppr[i] = prdmp[ic];
            elepr[i] = prele[ic];
        }
    }

    // process keywords containing atom specific Pauli repulsion
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "REPULSION") {
            ia = 0;
            spr = 0.;
            apr = 0.;
            epr = 0.;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> ia >> spr >> apr >> epr)) break;
            if (ia<0 and ia>=-n) {
                ia = -ia;
                ia--;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Pauli Repulsion Values for Specific Atoms :\n\n");
                    printf("        Atom                 Size            Damp        Valence\n\n");
                }
                if (!silent) {
                    printf("      %6d       %15.4f%15.4f%15.3f\n", ia+1, spr, apr, epr);
                }
                sizpr[ia] = spr;
                dmppr[ia] = apr;
                elepr[ia] = -std::abs(epr);
            }
        }
    }

    // condense repulsion sites to the list of multipole sites
    nrep = 0;
    if (use_repel) {
        for (int i = 0; i < n; i++) {
            if (sizpr[i] != 0.) {
                irep[nrep] = i;
                nrep++;
                replist[i] = nrep;
                for (int j = 0; j < maxpole; j++) {
                    repole[i][j] = pole[i][j];
                }
            }
        }
    }

    // test multipoles at chiral sites and invert if necessary
    chkrepole();

    // turn off the Pauli repulsion potential if not used
    if (nrep == 0) use_repel = false;
}
}
