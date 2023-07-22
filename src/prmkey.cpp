//////////////////////////////////////////////////////
//                                                  //
//  prmkey.cpp  --  interpret force field keywords  //
//                                                  //
//////////////////////////////////////////////////////

// "prmkey" parses a text string to extract keywords related to
// force field potential energy functional forms and constants


#include "angpot.h"
#include "bndpot.h"
#include "chgpot.h"
#include "ctrpot.h"
#include "dsppot.h"
#include "expol.h"
#include "extfld.h"
#include "fields.h"
#include "gettext.h"
#include "getword.h"
#include "mplpot.h"
#include "polpot.h"
#include "potent.h"
#include "prmkey.h"
#include "reppot.h"
#include "rxnpot.h"
#include "torpot.h"
#include "units.h"
#include "upcase.h"
#include "urypot.h"
#include "vdwpot.h"
#include <sstream>

void potoff();
void valoff();
void nbondoff();

void prmkey(std::string record)
{
    std::string value;
    std::string keyword;
    std::string string;
    std::istringstream iss;

    // parse the line to extract any possible keyword
    int next = 0;
    upcase(record);
    gettext(record,keyword,next);
    string = record.substr(next);
    iss.str(string);

    // select the individual force field potential terms
    if (keyword == "BONDTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_bond = true;
        if (value == "NONE") use_bond = false;
    }
    else if (keyword == "ANGLETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_angle = true;
        if (value == "NONE") use_angle = false;
    }
    else if (keyword == "STRBNDTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_strbnd = true;
        if (value == "NONE") use_strbnd = false;
    }
    else if (keyword == "UREYBRADTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_urey = true;
        if (value == "NONE") use_urey = false;
    }
    else if (keyword == "ANGANGTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_angang = true;
        if (value == "NONE") use_angang = false;
    }
    else if (keyword == "OPBENDTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_opbend = true;
        if (value == "NONE") use_opbend = false;
    }
    else if (keyword == "OPDISTTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_opdist = true;
        if (value == "NONE")  use_opdist = false;
    }
    else if (keyword == "IMPROPTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_improp = true;
        if (value == "NONE")  use_improp = false;
    }
    else if (keyword == "IMPTORTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_imptor = true;
        if (value == "NONE")  use_imptor = false;
    }
    else if (keyword == "TORSIONTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_tors = true;
        if (value == "NONE")  use_tors = false;
    }
    else if (keyword == "PITORSTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_pitors = true;
        if (value == "NONE")  use_pitors = false;
    }
    else if (keyword == "STRTORTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_strtor = true;
        if (value == "NONE")  use_strtor = false;
    }
    else if (keyword == "ANGTORTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_angtor = true;
        if (value == "NONE")  use_angtor = false;
    }
    else if (keyword == "TORTORTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_tortor = true;
        if (value == "NONE")  use_tortor = false;
    }
    else if (keyword == "VDWTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_vdw = true;
        if (value == "NONE")  use_vdw = false;
    }
    else if (keyword == "REPULSIONTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_repel = true;
        if (value == "NONE")  use_repel = false;
    }
    else if (keyword == "DISPERSIONTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_disp = true;
        if (value == "NONE")  use_disp = false;
    }
    else if (keyword == "CHARGETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_charge = true;
        if (value == "NONE")  use_charge = false;
    }
    else if (keyword == "CHGDPLTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_chgdpl = true;
        if (value == "NONE")  use_chgdpl = false;
    }
    else if (keyword == "DIPOLETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_dipole = true;
        if (value == "NONE")  use_dipole = false;
    }
    else if (keyword == "MULTIPOLETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_mpole = true;
        if (value == "NONE")  use_mpole = false;
    }
    else if (keyword == "POLARIZETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_polar = true;
        if (value == "NONE")  use_polar = false;
    }
    else if (keyword == "CHGTRNTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_chgtrn = true;
        if (value == "NONE")  use_chgtrn = false;
    }
    else if (keyword == "CHGFLXTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_chgflx = true;
        if (value == "NONE")  use_chgflx = false;
    }
    else if (keyword == "RXNFIELDTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_rxnfld = true;
        if (value == "NONE")  use_rxnfld = false;
    }
    else if (keyword == "SOLVATETERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_solv = true;
        if (value == "NONE")  use_solv = false;
    }
    else if (keyword == "METALTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_metal = true;
        if (value == "NONE")  use_metal = false;
    }
    else if (keyword == "RESTRAINTERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_geom = true;
        if (value == "NONE")  use_geom = false;
    }
    else if (keyword == "EXTRATERM") {
        getword(record,value,next);
        if (value == "ONLY") potoff();
        use_extra = true;
        if (value == "NONE")  use_extra = false;
    }
    else if (keyword == "VALENCETERM") {
        getword(record,value,next);
        if (value == "ONLY") nbondoff();
        if (value == "NONE") valoff();
    }
    else if (keyword == "NONBONDTERM") {
        getword(record,value,next);
        if (value == "ONLY") valoff();
        if (value == "NONE") nbondoff();
    }

    // select the name of the force field parameter set
    if (keyword == "FORCEFIELD") {
        getword (record,forcefield,next);
    }

    // set control parameters for bond stretching potentials
    else if (keyword == "BONDTYPE") {
        getword(record,bndtyp,next);
    }
    else if (keyword == "BONDUNIT") {
        iss >> bndunit;
    }
    else if (keyword == "BOND-CUBIC") {
        iss >> cbnd;
    }
    else if (keyword == "BOND-QUARTIC") {
        iss >> qbnd;
    }

    // set control parameters for bond angle bending potentials
    else if (keyword == "ANGLEUNIT") {
        iss >> angunit;
    }
    else if (keyword == "ANGLE-CUBIC") {
        iss >> cang;
    }
    else if (keyword == "ANGLE-QUARTIC") {
        iss >> qang;
    }
    else if (keyword == "ANGLE-PENTIC") {
        iss >> pang;
    }
    else if (keyword == "ANGLE-SEXTIC") {
        iss >> sang;
    }

    // set control parameters for stretch-bend potentials
    else if (keyword == "STRBNDUNIT") {
        iss >> stbnunit;
    }

    // set control parameters for Urey-Bradley potentials
    else if (keyword == "UREYUNIT") {
        iss >> ureyunit;
    }
    else if (keyword == "UREY-CUBIC") {
        iss >> cury;
    }
    else if (keyword == "UREY-QUARTIC") {
        iss >> qury;
    }

    // set control parameters for out-of-plane bend potentials
    else if (keyword == "OPBENDTYPE") {
        getword(record,opbtyp,next);
    }
    else if (keyword == "OPBENDUNIT") {
        iss >> opbunit;
    }
    else if (keyword == "OPBEND-CUBIC") {
        iss >> copb;
    }
    else if (keyword == "OPBEND-QUARTIC") {
        iss >> qopb;
    }
    else if (keyword == "OPBEND-PENTIC") {
        iss >> popb;
    }
    else if (keyword == "OPBEND-SEXTIC") {
        iss >> sopb;
    }

    // set control parameters for out-of-plane distance potentials
    else if (keyword == "OPDISTUNIT") {
        iss >> opdunit;
    }
    else if (keyword == "OPDIST-CUBIC") {
        iss >> copd;
    }
    else if (keyword == "OPDIST-QUARTIC") {
        iss >> qopd;
    }
    else if (keyword == "OPDIST-PENTIC") {
        iss >> popd;
    }
    else if (keyword == "OPDIST-SEXTIC") {
        iss >> sopd;
    }

    // set control parameters for other local geometry potentials
    else if (keyword == "ANGANGUNIT") {
        iss >> aaunit;
    }
    else if (keyword == "IMPROPUNIT") {
        iss >> idihunit;
    }
    else if (keyword == "IMPTORUNIT") {
        iss >> itorunit;
    }
    else if (keyword == "TORSIONUNIT") {
        iss >> torsunit;
    }
    else if (keyword == "PITORSUNIT") {
        iss >> ptorunit;
    }
    else if (keyword == "STRTORUNIT") {
        iss >> storunit;
    }
    else if (keyword == "ANGTORUNIT") {
        iss >> atorunit;
    }
    else if (keyword == "TORTORUNIT") {
        iss >> ttorunit;
    }

    // set control parameters for van der Waals potentials
    else if (keyword == "VDWINDEX") {
        getword(record,vdwindex,next);
    }
    else if (keyword == "VDWTYPE") {
        getword(record,vdwtyp,next);
    }
    else if (keyword == "RADIUSTYPE") {
        getword(record,radtyp,next);
    }
    else if (keyword == "RADIUSSIZE") {
        getword(record,radsiz,next);
    }
    else if (keyword == "RADIUSRULE") {
        getword(record,radrule,next);
    }
    else if (keyword == "EPSILONRULE") {
        getword(record,epsrule,next);
    }
    else if (keyword == "GAUSSTYPE") {
        getword(record,gausstyp,next);
    }
    else if (keyword == "A-EXPTERM") {
        iss >> abuck;
    }
    else if (keyword == "B-EXPTERM") {
        iss >> bbuck;
    }
    else if (keyword == "C-EXPTERM") {
        iss >> cbuck;
    }
    else if (keyword == "GAMMA-HALGREN") {
        iss >> ghal;
    }
    else if (keyword == "DELTA-HALGREN") {
        iss >> dhal;
    }
    else if (keyword == "VDW-12-SCALE") {
        iss >> v2scale;
        if (v2scale > 1.) v2scale = 1. / v2scale;
    }
    else if (keyword == "VDW-13-SCALE") {
        iss >> v3scale;
        if (v3scale > 1.) v3scale = 1. / v3scale;
    }
    else if (keyword == "VDW-14-SCALE") {
        iss >> v4scale;
        if (v4scale > 1.) v4scale = 1. / v4scale;
    }
    else if (keyword == "VDW-15-SCALE") {
        iss >> v5scale;
        if (v5scale > 1.) v5scale = 1. / v5scale;
    }
    else if (keyword == "VDW-CORRECTION") {
        use_vcorr = true;
    }

    // set control parameters for Pauli repulsion potential
    else if (keyword == "REP-12-SCALE") {
        iss >> r2scale;
        if (r2scale > 1.) r2scale = 1. / r2scale;
    }
    else if (keyword == "REP-13-SCALE") {
        iss >> r3scale;
        if (r3scale > 1.) r3scale = 1. / r3scale;
    }
    else if (keyword == "REP-14-SCALE") {
        iss >> r4scale;
        if (r4scale > 1.) r4scale = 1. / r4scale;
    }
    else if (keyword == "REP-15-SCALE") {
        iss >> r5scale;
        if (r5scale > 1.) r5scale = 1. / r5scale;
    }

    // set control parameters for dispersion potential
    else if (keyword == "DISP-12-SCALE") {
        iss >> dsp2scale;
        if (dsp2scale > 1.) dsp2scale = 1. / dsp2scale;
    }
    else if (keyword == "DISP-13-SCALE") {
        iss >> dsp3scale;
        if (dsp3scale > 1.) dsp3scale = 1. / dsp3scale;
    }
    else if (keyword == "DISP-14-SCALE") {
        iss >> dsp4scale;
        if (dsp4scale > 1.) dsp4scale = 1. / dsp4scale;
    }
    else if (keyword == "DISP-15-SCALE") {
        iss >> dsp5scale;
        if (dsp5scale > 1.) dsp5scale = 1. / dsp5scale;
    }
    else if (keyword == "DISP-CORRECTION") {
        use_dcorr = true;
    }

    // set control parameters for charge-charge potentials
    else if (keyword == "ELECTRIC") {
        iss >> electric;
    }
    else if (keyword == "DIELECTRIC") {
        iss >> dielec;
    }
    else if (keyword == "CHG-BUFFER") {
        iss >> ebuffer;
    }
    else if (keyword == "CHG-11-SCALE") {
        iss >> c1scale;
        if (c1scale > 1.) c1scale = 1. / c1scale;
    }
    else if (keyword == "CHG-12-SCALE") {
        iss >> c2scale;
        if (c2scale > 1.) c2scale = 1. / c2scale;
    }
    else if (keyword == "CHG-13-SCALE") {
        iss >> c3scale;
        if (c3scale > 1.) c3scale = 1. / c3scale;
    }
    else if (keyword == "CHG-14-SCALE") {
        iss >> c4scale;
        if (c4scale > 1.) c4scale = 1. / c4scale;
    }
    else if (keyword == "CHG-15-SCALE") {
        iss >> c5scale;
        if (c5scale > 1.) c5scale = 1. / c5scale;
    }
    else if (keyword == "NEIGHBOR-GROUPS") {
        neutnbr = true;
    }
    else if (keyword == "NEUTRAL-GROUPS") {
        neutcut = true;
    }
    else if (keyword == "EXTERNAL-FIELD") {
        iss >> exfld[0] >> exfld[1] >> exfld[2];
        use_exfld = true;
        for (int i = 0; i < 3; i++) {
            exfld[i] = exfld[i] / elefield;
        }
    }

    // set control parameters for atomic multipole potentials
    else if (keyword == "PENETRATION") {
        getword(record,pentyp,next);
    }
    else if (keyword == "MPOLE-12-SCALE") {
        iss >> m2scale;
        if (m2scale > 1.) m2scale = 1. / m2scale;
    }
    else if (keyword == "MPOLE-13-SCALE") {
        iss >> m3scale;
        if (m3scale > 1.) m3scale = 1. / m3scale;
    }
    else if (keyword == "MPOLE-14-SCALE") {
        iss >> m4scale;
        if (m4scale > 1.) m4scale = 1. / m4scale;
    }
    else if (keyword == "MPOLE-15-SCALE") {
        iss >> m5scale;
        if (m5scale > 1.) m5scale = 1. / m5scale;
    }

    // set control parameters for polarization potentials
    else if (keyword == "POLARIZATION") {
        getword(record,poltyp,next);
    }
    else if (keyword == "EXCHANGE-POLAR") {
        getword(record,scrtyp,next);
    }
    else if (keyword == "POLAR-ITER") {
        iss >> politer;
    }
    else if (keyword == "POLAR-EPS") {
        iss >> poleps;
    }
    else if (keyword == "USOLVE-ACCEL") {
        iss >> uaccel;
    }
    else if (keyword == "D-EQUALS-P") {
        dpequal = true;
    }
    else if (keyword == "POLAR-12-SCALE") {
        iss >> p2scale;
        if (p2scale > 1.) p2scale = 1. / p2scale;
    }
    else if (keyword == "POLAR-13-SCALE") {
        iss >> p3scale;
        if (p3scale > 1.) p3scale = 1. / p3scale;
    }
    else if (keyword == "POLAR-14-SCALE") {
        iss >> p4scale;
        if (p4scale > 1.) p4scale = 1. / p4scale;
    }
    else if (keyword == "POLAR-15-SCALE") {
        iss >> p5scale;
        if (p5scale > 1.) p5scale = 1. / p5scale;
    }
    else if (keyword == "POLAR-12-INTRA") {
        iss >> p2iscale;
        if (p2iscale > 1.) p2iscale = 1. / p2iscale;
    }
    else if (keyword == "POLAR-13-INTRA") {
        iss >> p3iscale;
        if (p3iscale > 1.) p3iscale = 1. / p3iscale;
    }
    else if (keyword == "POLAR-14-INTRA") {
        iss >> p4iscale;
        if (p4iscale > 1.) p4iscale = 1. / p4iscale;
    }
    else if (keyword == "POLAR-15-INTRA") {
        iss >> p5iscale;
        if (p5iscale > 1.) p5iscale = 1. / p5iscale;
    }
    else if (keyword == "DIRECT-11-SCALE") {
        iss >> d1scale;
        if (d1scale > 1.) d1scale = 1. / d1scale;
    }
    else if (keyword == "DIRECT-12-SCALE") {
        iss >> d2scale;
        if (d2scale > 1.) d2scale = 1. / d2scale;
    }
    else if (keyword == "DIRECT-13-SCALE") {
        iss >> d3scale;
        if (d3scale > 1.) d3scale = 1. / d3scale;
    }
    else if (keyword == "DIRECT-14-SCALE") {
        iss >> d4scale;
        if (d4scale > 1.) d4scale = 1. / d4scale;
    }
    else if (keyword == "MUTUAL-11-SCALE") {
        iss >> u1scale;
        if (u1scale > 1.) u1scale = 1. / u1scale;
    }
    else if (keyword == "MUTUAL-12-SCALE") {
        iss >> u2scale;
        if (u2scale > 1.) u2scale = 1. / u2scale;
    }
    else if (keyword == "MUTUAL-13-SCALE") {
        iss >> u3scale;
        if (u3scale > 1.) u3scale = 1. / u3scale;
    }
    else if (keyword == "MUTUAL-14-SCALE") {
        iss >> u4scale;
        if (u4scale > 1.) u4scale = 1. / u4scale;
    }
    else if (keyword == "INDUCE-12-SCALE") {
        iss >> w2scale;
        if (w2scale > 1.) w2scale = 1. / w2scale;
    }
    else if (keyword == "INDUCE-13-SCALE") {
        iss >> w3scale;
        if (w3scale > 1.) w3scale = 1. / w3scale;
    }
    else if (keyword == "INDUCE-14-SCALE") {
        iss >> w4scale;
        if (w4scale > 1.) w4scale = 1. / w4scale;
    }
    else if (keyword == "INDUCE-15-SCALE") {
        iss >> w5scale;
        if (w5scale > 1.) w5scale = 1. / w5scale;
    }

    // set control parameters for charge transfer potentials
    else if (keyword == "CHARGETRANSFER ") {
        getword(record,ctrntyp,next);
    }

    // set control parameters for reaction field potentials
    else if (keyword == "REACTIONFIELD ") {
        iss >> rfsize >> rfbulkd >> rfterms;
    }
}


////////////////////////////////////////////////////////
//                                                    //
//  potoff.cpp  --  turn off all potential functions  //
//                                                    //
////////////////////////////////////////////////////////

// "potoff" clears the forcefield definition by turning off
// the use of each of the potential energy functions


void potoff()
{
    // turn off the use of each of the potential energy functions
    use_bond = false;
    use_angle = false;
    use_strbnd = false;
    use_urey = false;
    use_angang = false;
    use_opbend = false;
    use_opdist = false;
    use_improp = false;
    use_imptor = false;
    use_tors = false;
    use_pitors = false;
    use_strtor = false;
    use_angtor = false;
    use_tortor = false;
    use_vdw = false;
    use_repel = false;
    use_disp = false;
    use_charge = false;
    use_chgdpl = false;
    use_dipole = false;
    use_mpole = false;
    use_polar = false;
    use_chgtrn = false;
    use_rxnfld = false;
    use_solv = false;
    use_metal = false;
    use_geom = false;
    use_extra = false;
}


////////////////////////////////////////////////////////
//                                                    //
//  valoff.cpp  --  turn off valence potential terms  //
//                                                    //
////////////////////////////////////////////////////////

// "valoff" turns off the use of each of the valence
// potential energy functions


void valoff()
{
    // turn off the use of each of the valence energy functions
    use_bond = false;
    use_angle = false;
    use_strbnd = false;
    use_urey = false;
    use_angang = false;
    use_opbend = false;
    use_opdist = false;
    use_improp = false;
    use_imptor = false;
    use_tors = false;
    use_pitors = false;
    use_strtor = false;
    use_angtor = false;
    use_tortor = false;
    use_geom = false;
}


//////////////////////////////////////////////////////////
//                                                      //
//  nbondoff.cpp  --  turn off nonbond potential terms  //
//                                                      //
//////////////////////////////////////////////////////////

// "nbondoff" turns off the use of each of the nonbonded
// potential energy functions


void nbondoff()
{
    // turn off the use of each of the nonbonded energy functions
    use_vdw = false;
    use_repel = false;
    use_disp = false;
    use_charge = false;
    use_chgdpl = false;
    use_dipole = false;
    use_mpole = false;
    use_polar = false;
    use_chgtrn = false;
    use_rxnfld = false;
    use_solv = false;
    use_metal = false;
}
