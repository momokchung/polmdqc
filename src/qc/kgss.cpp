// Author: Moses KJ Chung
// Year:   2024

#include "gettext.h"
#include "gss.h"
#include "keys.h"
#include "kgss.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  kgss  --  guess parameter assignment  //
//                                        //
////////////////////////////////////////////

// "kgss" assigns guess parameters for quantum calculations

void kgss()
{
    int next;
    std::string guess;
    std::string keyword,record,string;
    std::istringstream iss;

    // process keywords containing guess parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "BASIS-GUESS") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            guess = "";
            iss >> guess;
            if (guess == "TRUE") bssguess = true;
            else if (guess == "FALSE") bssguess = false;
        }
        else if (keyword == "DENSITY-GUESS") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            guess = "";
            iss >> guess;
            if (guess == "SAD") denguess = DenGuess::SAD;
            else if (guess == "CORE") denguess = DenGuess::Core;
            else if (guess == "SAP") denguess = DenGuess::SAP;
        }
    }
}
}
