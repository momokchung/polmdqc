/////////////////////////////////////////////////////////
//                                                     //
//  field.cpp  --  get the potential energy functions  //
//                                                     //
/////////////////////////////////////////////////////////

// "field" sets the force field potential energy functions from
// a parameter file and modifications specified in a keyfile


#include "fatal.h"
#include "field.h"
#include "fields.h"
#include "getprm.h"
#include "getstring.h"
#include "gettext.h"
#include "getword.h"
#include "inform.h"
#include "keys.h"
#include "potent.h"
#include "prmkey.h"
#include "sizes.h"
#include "upcase.h"
#include <sstream>

void field()
{
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;
    int ia,ib,next;
    bool header;

    // set the default values for the active potentials
    use_bond = true;
    use_angle = true;
    use_strbnd = true;
    use_urey = true;
    use_angang = true;
    use_opbend = true;
    use_opdist = true;
    use_improp = true;
    use_imptor = true;
    use_tors = true;
    use_pitors = true;
    use_strtor = true;
    use_angtor = true;
    use_tortor = true;
    use_vdw = true;
    use_repel = true;
    use_disp = true;
    use_charge = true;
    use_chgdpl = true;
    use_dipole = true;
    use_mpole = true;
    use_polar = true;
    use_chgtrn = true;
    use_chgflx = true;
    use_rxnfld = false;
    use_solv = true;
    use_metal = false;
    use_geom = true;
    use_extra = true;

    // read the potential energy force field parameter file
    getprm();

    // check keywords for biopolymer atom type definitions
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "BIOTYPE") {
            ia = 0;
            ib = 0;
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ia;
            getword(record,string,next);
            getstring(record,string,next);
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ib;
            if (ia>=0 and ia<=maxbio) {
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Biopolymer Type Definitions\n\n     BioType          Atom Type\n\n");
                }
                biotyp[ia-1] = ib-1;
                if (!silent) {
                    printf(" %8d        %8d\n", ia,ib);
                }
            }
            else if (ia > maxbio) {
                printf("\n FIELDS  --  Too many Biopolymer Types; Increase MAXBIO\n");
                fatal();
            }
        }
    }

    // check keywords for potential function control parameters
    for (int i = 0; i < nkey; i++) {
        record = keyline[i];
        prmkey(record);
    }
}
