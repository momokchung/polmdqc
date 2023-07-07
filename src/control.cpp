/////////////////////////////////////////////////////////
//                                                     //
//  control.cpp  --  set information and output types  //
//                                                     //
/////////////////////////////////////////////////////////

// "control" gets initial values for parameters that determine
// the output style and information level provided by PolQCMD


#include "argue.h"
#include "control.h"
#include "gettext.h"
#include "inform.h"
#include "keys.h"
#include "output.h"
#include "upcase.h"

void control()
{
    // set default values for information and output variables
    digits = 4;
    verbose = false;
    debug = false;
    holdup = false;
    informAbort = false;
    arcsave = false;
    dcdsave = false;
    cyclesave = false;
    noversion = false;
    overwrite = false;

    // check for control parameters on the command line
    bool exist = false;
    for (int i = 0; i < narg-1; i++) {
        std::string string = arg[i];
        upcase(string);
        if (string.substr(0,2).compare("-D") == 0) {
            debug = true;
            verbose = true;
        }
        else if (string.substr(0,2).compare("-V") == 0) {
        verbose = true;
        }
    }

    // search keywords for various control parameters
    for (int i = 0; i < nkey; i++) {
        int next = 0;
        std::string record = keyline[i];
        std::string keyword;
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "DIGITS") {
            gettext(record,keyword,next);
            try {
                digits = std::stoi(keyword);
            } catch (const std::exception& e) {}
        }
        else if (keyword == "DEBUG") {
            debug = true;
            verbose = true;
        }
        else if (keyword == "VERBOSE") {
            verbose = true;
        }
        else if (keyword == "EXIT-PAUSE") {
            holdup = true;
        }
        else if (keyword == "ARCHIVE") {
            arcsave = true;
            dcdsave = false;
        }
        else if (keyword == "DCD-ARCHIVE") {
            arcsave = false;
            dcdsave = true;
        }
        else if (keyword == "NOARCHIVE") {
            arcsave = false;
            dcdsave = false;
            cyclesave = false;
        }
        else if (keyword == "SAVE-CYCLE") {
            dcdsave = false;
            cyclesave = true;
        }
        else if (keyword == "NOVERSION") {
            noversion = true;
        }
        else if (keyword == "OVERWRITE") {
            overwrite = true;
        }
    }
}