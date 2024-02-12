// Author: Moses KJ Chung
// Year:   2024

#include "gettext.h"
#include "keys.h"
#include "kscf.h"
#include "scft.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////
//                                   //
//  kscf  --  assign scf parameters  //
//                                   //
///////////////////////////////////////

// "kscf" assigns scf parameters for quantum calculations

void kscf()
{
    int next;
    std::string scfs;
    std::string keyword,record,string;
    std::istringstream iss;

    // process keywords containing scf parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "SCF-TYPE") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            scfs = "";
            iss >> scfs;
            if (scfs == "PK") scftyp = SCFType::PK;
            else if (scfs == "DF") scftyp = SCFType::DF;
            else if (scfs == "DIRECT") scftyp = SCFType::Direct;
        }
    }
}
}
