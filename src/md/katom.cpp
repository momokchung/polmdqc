// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "couple.h"
#include "getnumb.h"
#include "getstring.h"
#include "gettext.h"
#include "inform.h"
#include "katom.h"
#include "katoms.h"
#include "keys.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  katom  --  atom type parameter assignment  //
//                                             //
/////////////////////////////////////////////////

// "katom" assigns an atom type definitions to each atom in
// the structure and processes any new or changed values

// literature reference:

// K. A. Feenstra, B. Hess and H. J. C. Berendsen, "Improving
// Efficiency of Large Time-Scale Molecular Dynamics Simulations
// of Hydrogen-Rich Systems", Journal of Computational Chemistry,
// 8, 786-798 (1999)

// C. W. Hopkins, S. Le Grand, R. C. Walker and A. E. Roitberg,
// "Long-Time-Step Molecular Dynamics through Hydrogen Mass
// Repartitioning", Journal of Chemical Theory and Computation,
// 11, 1864-1874 (2015)

void katom()
{
    int next,nh;
    int cls,atn,lig;
    double wght,sum;
    double hmax,hmass;
    double dmin,dmass;
    bool header,heavy;
    std::string symb,notice;
    std::string keyword,record,string;
    std::istringstream iss;

    // process keywords containing atom type parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ATOM") {
            int k = 0;
            cls = 0;
            symb = "";
            notice = "";
            atn = 0;
            wght = 0.;
            lig = 0;
            getnumb(record,k,next);
            if (k>0 and k<=maxtyp) {
                getnumb(record,cls,next);
                if (cls == 0)  cls = k;
                int km = k-1;
                int clsm = cls-1;
                atmcls[km] = clsm;
                gettext(record,symb,next);
                getstring(record,notice,next);
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> atn >> wght >> lig;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Atom Definition Parameters :\n\n");
                    printf("     Type  Class  Symbol  Description");
                    printf("               Atomic    Mass   Valence\n\n");
                }
                symbol[km] = symb;
                describe[km] = notice;
                atmnum[km] = atn;
                weight[km] = wght;
                ligand[km] = lig;
                if (!silent) {
                    printf(" %8d%6d     %-3s   %-24s%6d%11.3f%6d\n", k, cls, symb.c_str(), notice.c_str(), atn, wght, lig);
                }
            }
            else if (k > maxtyp) {
                printf("\n KATOM   --  Too many Atom Types; Increase MAXTYP\n");
                informAbort = true;
            }
        }
    }

    // transfer atom type values to individual atoms
    for (int i = 0; i < n; i++) {
        int k = type[i];
        if (k == -1 or k >= maxtyp) {
            atomClass[i] = -1;
            atomic[i] = 0;
            mass[i] = 0.;
            valence[i] = 0;
            story[i] = "Undefined Atom Type";
        }
        else {
            if (symbol[k] != "")  name[i] = symbol[k];
            atomClass[i] = atmcls[k];
            atomic[i] = atmnum[k];
            mass[i] = weight[k];
            valence[i] = ligand[k];
            story[i] = describe[k];
        }
    }

    // repartition hydrogen masses to use "heavy" hydrogens
    heavy = false;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "HEAVY-HYDROGEN") {
            heavy = true;
            hmax = 4.;
            iss >> hmax;
        }
    }
    if (heavy) {
        for (int i = 0; i < n; i++) {
            nh = 0;
            sum = mass[i];
            for (int j = 0; j < n12[i]; j++) {
                int k = i12[i][j];
                if (atomic[k] == 1) {
                    nh++;
                    sum += mass[k];
                }
            }
            hmass = std::min(hmax,sum/static_cast<double>(nh+1));
            for (int j = 0; j < n12[i]; j++) {
                int k = i12[i][j];
                if (atomic[k] == 1) {
                    dmass = hmass - mass[k];
                    mass[k] = hmass;
                    mass[i] -= dmass;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            if (mass[i] < hmax) {
                dmass = hmax - mass[i];
                dmin = hmax + dmass;
                for (int j = 0; j < n12[i]; j++) {
                    int k = i12[i][j];
                    if (mass[k] > dmin) {
                        mass[k] -= dmass;
                        mass[i] = hmax;
                        goto label_60;
                    }
                }
                for (int j = 0; j < n13[i]; j++) {
                    int k = i13[i][j];
                    if (mass[k] > dmin) {
                        mass[k] -= dmass;
                        mass[i] = hmax;
                        goto label_60;
                    }
                }
                label_60:
                continue;
            }
        }
    }

    // process keywords containing atom types for specific atoms
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ATOM") {
            int k = 0;
            symb = "";
            notice = "";
            atn = 0;
            wght = 0.;
            lig = 0;
            getnumb(record,k,next);
            if (k<0 and k>=-n) {
                getnumb(record,cls,next);
                gettext(record,symb,next);
                getstring(record,notice,next);
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> atn >> wght >> lig;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Atom Definitions for Specific Atoms :\n\n");
                    printf("     Atom  Class  Symbol  Description");
                    printf("               Atomic    Mass   Valence\n\n");
                }
                k = -k;
                if (cls == 0)  cls = k;
                int km = k-1;
                int clsm = cls-1;
                atomClass[k] = clsm;
                name[k] = symb;
                story[k] = notice;
                atomic[k] = atn;
                mass[k] = wght;
                valence[k] = lig;
                if (!silent) {
                    printf(" %8d%6d     %-3s   %-24s%6d%11.3f%6d\n", k, cls, symb.c_str(), notice.c_str(), atn, wght, lig);
                }
            }
        }
    }

    // check for presence of undefined atom types or classes
    header = true;
    for (int i = 0; i < n; i++) {
        int k = type[i];
        cls = atomClass[i];
        if (k<0 or k>maxtyp-1 or cls<0 or cls>maxclass-1) {
            informAbort = true;
            if (header) {
                header = false;
                printf("\n Undefined Atom Types or Classes :\n\n");
                printf(" Type          Atom Number     Atom Type     Atom Class\n\n");
            }
            printf(" Atom         %8d          %5d          %5d\n", i+1, k+1, cls+1);
        }
    }

    // check the number of atoms attached to each atom
    header = true;
    for (int i = 0; i < n; i++) {
        if (n12[i] != valence[i]) {
            if (header) {
                header = false;
                printf("\n Atoms with an Unusual Number of Attached Atoms :\n\n");
                printf(" Type           Atom Name      Atom Type       Expected    Found\n\n");
            }
            printf(" Valence    %8d-%-3s        %5d          %5d     %5d\n", i+1, name[i].c_str(), type[i]+1, valence[i], n12[i]);
        }
    }
}
}
