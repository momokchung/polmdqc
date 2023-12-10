// Author: Moses KJ Chung
// Year:   2023

#include "action.h"
#include "energi.h"
#include "inform.h"
#include "mdqclimits.h"
#include "partyze.h"
#include "potent.h"
#include <cmath>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  partyze  --  energy component decomposition  //
//                                               //
///////////////////////////////////////////////////

// "partyze" prints the energy component and number of
// interactions for each of the potential energy terms

void partyze()
{
    int numSpaces;
    int width;
    int precision;
    int intWidth;
    int tempSpaces;

    // write out each energy component to the desired precision
    numSpaces = 5;
    width = 16;
    precision = 4;
    intWidth = 17;
    if (digits >= 6) {
        numSpaces = 3;
        width = 18;
        precision = 6;
    }
    if (digits >= 8) {
        numSpaces = 1;
        width = 20;
        precision = 8;
    }
    printf("\n Energy Component Breakdown :           Kcal/mole        Interactions\n\n");
    if (use_bond and (neb!=0 or eb!=0.)) {
        tempSpaces = numSpaces + 12;
        printf(" Bond Stretching%*s%*.*f%*d\n", tempSpaces, "", width, precision, eb, intWidth, neb);
    }
    if (use_angle and (nea!=0 or ea!=0.)) {
        tempSpaces = numSpaces + 14;
        printf(" Angle Bending%*s%*.*f%*d\n", tempSpaces, "", width, precision, ea, intWidth, nea);
    }
    if (use_strbnd and (neba!=0 or eba!=0.)) {
        tempSpaces = numSpaces + 15;
        printf(" Stretch-Bend%*s%*.*f%*d\n", tempSpaces, "", width, precision, eba, intWidth, neba);
    }
    if (use_urey and (neub!=0 or eub!=0.)) {
        tempSpaces = numSpaces + 15;
        printf(" Urey-Bradley%*s%*.*f%*d\n", tempSpaces, "", width, precision, eub, intWidth, neub);
    }
    if (use_angang and (neaa!=0 or eaa!=0.)) {
        tempSpaces = numSpaces + 16;
        printf(" Angle-Angle%*s%*.*f%*d\n", tempSpaces, "", width, precision, eaa, intWidth, neaa);
    }
    if (use_opbend and (neopb!=0 or eopb!=0.)) {
        tempSpaces = numSpaces + 10;
        printf(" Out-of-Plane Bend%*s%*.*f%*d\n", tempSpaces, "", width, precision, eopb, intWidth, neopb);
    }
    if (use_opdist and (neopd!=0 or eopd!=0.)) {
        tempSpaces = numSpaces + 6;
        printf(" Out-of-Plane Distance%*s%*.*f%*d\n", tempSpaces, "", width, precision, eopd, intWidth, neopd);
    }
    if (use_improp and (neid!=0 or eid!=0.)) {
        tempSpaces = numSpaces + 10;
        printf(" Improper Dihedral%*s%*.*f%*d\n", tempSpaces, "", width, precision, eid, intWidth, neid);
    }
    if (use_imptor and (neit!=0 or eit!=0.)) {
        tempSpaces = numSpaces + 11;
        printf(" Improper Torsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, eit, intWidth, neit);
    }
    if (use_tors and (net!=0 or et!=0.)) {
        tempSpaces = numSpaces + 12;
        printf(" Torsional Angle%*s%*.*f%*d\n", tempSpaces, "", width, precision, et, intWidth, net);
    }
    if (use_pitors and (nept!=0 or ept!=0.)) {
        tempSpaces = numSpaces + 9;
        printf(" Pi-Orbital Torsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, ept, intWidth, nept);
    }
    if (use_strtor and (nebt!=0 or ebt!=0.)) {
        tempSpaces = numSpaces + 12;
        printf(" Stretch-Torsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, ebt, intWidth, nebt);
    }
    if (use_angtor and (neat!=0 or eat!=0.)) {
        tempSpaces = numSpaces + 14;
        printf(" Angle-Torsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, eat, intWidth, neat);
    }
    if (use_tortor and (nett!=0 or ett!=0.)) {
        tempSpaces = numSpaces + 12;
        printf(" Torsion-Torsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, ett, intWidth, nett);
    }

    if (use_vdw and (nev!=0 or ev!=0.)) {
        tempSpaces = numSpaces + 14;
        if (std::abs(ev) < 1e10) {
            printf(" Van der Waals%*s%*.*f%*d\n", tempSpaces, "", width, precision, ev, intWidth, nev);
        }
        else {
            printf(" Van der Waals%*s%*.*e%*d\n", tempSpaces, "", width, precision, ev, intWidth, nev);
        }
    }
    if (use_repel and (ner!=0 or er!=0.)) {
        tempSpaces = numSpaces + 18;
        if (std::abs(er) < 1e10) {
            printf(" Repulsion%*s%*.*f%*d\n", tempSpaces, "", width, precision, er, intWidth, ner);
        }
        else {
            printf(" Repulsion%*s%*.*e%*d\n", tempSpaces, "", width, precision, er, intWidth, ner);
        }
    }
    if (use_disp and (nedsp!=0 or edsp!=0.)) {
        tempSpaces = numSpaces + 17;
        printf(" Dispersion%*s%*.*f%*d\n", tempSpaces, "", width, precision, edsp, intWidth, nedsp);
    }
    if (use_charge and (nec!=0 or ec!=0.)) {
        tempSpaces = numSpaces + 14;
        if (std::abs(ec) < 1e10) {
            printf(" Charge-Charge%*s%*.*f%*d\n", tempSpaces, "", width, precision, ec, intWidth, nec);
        }
        else {
            printf(" Charge-Charge%*s%*.*e%*d\n", tempSpaces, "", width, precision, ec, intWidth, nec);
        }
    }
    if (use_chgdpl and (necd!=0 or ecd!=0.)) {
        tempSpaces = numSpaces + 14;
        if (std::abs(ecd) < 1e10) {
            printf(" Charge-Dipole%*s%*.*f%*d\n", tempSpaces, "", width, precision, ecd, intWidth, necd);
        }
        else {
            printf(" Charge-Dipole%*s%*.*e%*d\n", tempSpaces, "", width, precision, ecd, intWidth, necd);
        }
        
    }
    if (use_dipole and (ned!=0 or ed!=0.)) {
        tempSpaces = numSpaces + 14;
        if (std::abs(ed) < 1e10) {
            printf(" Dipole-Dipole%*s%*.*f%*d\n", tempSpaces, "", width, precision, ed, intWidth, ned);
        }
        else {
            printf(" Dipole-Dipole%*s%*.*e%*d\n", tempSpaces, "", width, precision, ed, intWidth, ned);
        }
        
    }
    if (use_mpole and (nem!=0 or em!=0.)) {
        tempSpaces = numSpaces + 10;
        if (std::abs(em) < 1e10) {
            printf(" Atomic Multipoles%*s%*.*f%*d\n", tempSpaces, "", width, precision, em, intWidth, nem);
        }
        else {
            printf(" Atomic Multipoles%*s%*.*e%*d\n", tempSpaces, "", width, precision, em, intWidth, nem);
        }
        
    }
    if (use_polar and (nep!=0 or ep!=0.)) {
        tempSpaces = numSpaces + 15;
        if (std::abs(ep) < 1e10) {
            printf(" Polarization%*s%*.*f%*d\n", tempSpaces, "", width, precision, ep, intWidth, nep);
        }
        else {
            printf(" Polarization%*s%*.*e%*d\n", tempSpaces, "", width, precision, ep, intWidth, nep);
        }
    }
    if (use_chgtrn and (nect!=0 or ect!=0.)) {
        tempSpaces = numSpaces + 12;
        if (std::abs(ect) < 1e10) {
            printf(" Charge Transfer%*s%*.*f%*d\n", tempSpaces, "", width, precision, ect, intWidth, nect);
        }
        else {
            printf(" Charge Transfer%*s%*.*e%*d\n", tempSpaces, "", width, precision, ect, intWidth, nect);
        }
        
    }
    if (use_rxnfld and (nerxf!=0 or erxf!=0.)) {
        tempSpaces = numSpaces + 13;
        printf(" Reaction Field%*s%*.*f%*d\n", tempSpaces, "", width, precision, erxf, intWidth, nerxf);
    }
    if (use_solv and (nes!=0 or es!=0.)) {
        tempSpaces = numSpaces + 9;
        printf(" Implicit Solvation%*s%*.*f%*d\n", tempSpaces, "", width, precision, es, intWidth, nes);
    }
    if (use_metal and (nelf!=0 or elf!=0.)) {
        tempSpaces = numSpaces + 9;
        printf(" Metal Ligand Field%*s%*.*f%*d\n", tempSpaces, "", width, precision, elf, intWidth, nelf);
    }
    if (use_geom and (neg!=0 or eg!=0.)) {
        tempSpaces = numSpaces + 7;
        printf(" Geometric Restraints%*s%*.*f%*d\n", tempSpaces, "", width, precision, eg, intWidth, neg);
    }
    if (use_extra and (nex!=0 or ex!=0.)) {
        tempSpaces = numSpaces + 9;
        printf(" Extra Energy Terms%*s%*.*f%*d\n", tempSpaces, "", width, precision, ex, intWidth, nex);
    }
}
}
