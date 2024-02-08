// Author: Moses KJ Chung
// Year:   2023

#include "inform.h"
#include "gettext.h"
#include "keys.h"
#include "optinit.h"
#include "output.h"
#include "potent.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  optinit  --  initialize structure optimization  //
//                                                  //
//////////////////////////////////////////////////////

// "optinit" initializes values and keywords used by multiple
// structure optimization methods

void optinit()
{
    int i,next;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // set default values for optimization parameters
    iprint = 0;
    iwrite = 0;
    frcsave = false;
    uindsave = false;

    // check for keywords containing any altered parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "PRINTOUT") {
            iss >> iprint;
        }
        else if (keyword == "WRITEOUT") {
            iss >> iwrite;
        }
        else if (keyword == "SAVE-FORCE") {
            frcsave = true;
        }
        else if (keyword == "SAVE-INDUCED") {
            uindsave = true;
        }
    }

    // check for use of induced dipole prediction methods
    // if (use_polar) predict();
}
}
