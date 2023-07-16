//////////////////////////////////////////////////////
//                                                  //
//  getprm.cpp  --  get force field parameter file  //
//                                                  //
//////////////////////////////////////////////////////

// "getprm" finds the potential energy parameter file
// and then opens and reads the parameters


#include "fatal.h"
#include "files.h"
#include "gettext.h"
#include "getprm.h"
#include "inform.h"
#include "initprm.h"
#include "inquire.h"
#include "keys.h"
#include "nextarg.h"
#include "params.h"
#include "setprm.h"
#include "suffix.h"
#include "upcase.h"
#include <iostream>
#include <fstream>
#include <sstream>

void getprm()
{
    int nask,next;
    bool exist,useprm;
    std::string  none;
    std::string  keyword;
    std::string  prmfile;
    std::string  prefix;
    std::string  record;
    std::string  string;
    std::istringstream iss;

    // set the default name for the parameter file
    useprm = true;
    prmfile = filename.substr(0,leng) + ".prm";

    // search the keyword list for the parameter filename
    for (int i = 0; i < nkey; i++) {
        next = 1;
        record = keyline[i];
        iss.clear();
        keyword = "";
        iss.str(record);
        iss >> keyword;
        upcase(keyword);
        if (keyword == "PARAMETERS" or keyword == "PARAMETER") {
            iss >> prmfile;
        }
    }
    // account for home directory abbreviation in filename
    if (prmfile.substr(0,2) == "~/") {
        char* homeDir = getenv("HOME");
        if (homeDir != nullptr) {
            std::string prefix = std::string(homeDir);
            prmfile = prefix + prmfile.substr(1);
        }
    }

    // check existence of default or specified parameter file
    suffix(prmfile, "prm", "old");
    exist = inquire(prmfile);

    // test for user specified absence of a parameter file
    if (!exist) {
        none = prmfile.substr(0,4);
        upcase(none);
        if (none == "NONE") {
            exist = true;
            useprm = false;
        }
    }

    // try to get a parameter filename from the command line
    if (!exist) {
        nextarg (prmfile, exist);
        if (exist) {
            suffix(prmfile, "prm", "old");
            exist = inquire(prmfile);
        }
    }

    // if necessary, ask for the parameter filename
    nask = 0;
    while (!exist and nask<maxask) {
        nask++;
        printf("\n Enter Parameter File Name [<Enter>=NONE] :  ");
        std::getline(std::cin, prmfile);
        iss.clear();
        none = "";
        iss.str(prmfile);
        bool input;
        if (!(iss >> none)) {
            input = false;
        }
        else {
            input = true;
            upcase(none);
        }
        if (!input) {
            exist = true;
            useprm = false;
        }
        else if (none == "NONE") {
            exist = true;
            useprm = false;
        }
        else {
            if (prmfile.substr(0,2) == "~/") {
                char* homeDir = getenv("HOME");
                if (homeDir != nullptr) {
                    std::string prefix = std::string(homeDir);
                    prmfile = prefix + prmfile.substr(1);
                }
            }
            suffix(prmfile, "prm", "old");
            exist = inquire(prmfile);
        }
    }
    if (!exist) fatal();

    // read the parameter file and store it for latter use
    nprm = 0;
    if (useprm) {
        std::ifstream file(prmfile);
        while (std::getline(file, record)) {
            prmline[nprm] = record;
            nprm++;
            if (nprm >= maxprm) {
                printf("\n GETPRM  --  Parameter File Too Large; Increase MAXPRM\n");
                fatal();
            }
        }
        file.close();
    }

    // convert underbar characters to dashes in all keywords
    for (int i = 0; i < nprm; i++) {
        int next = 0;
        record = prmline[i];
        gettext(record, keyword, next);
        for (int j = 0; j < next; j++)
        {
            if (record[j] == '_') record[j] = '-';
        }
        prmline[i] = record;
    }

    // count and allocate memory for the parameter values
    setprm();

    // initialize force field control and parameter values
    initprm();

    // get control and parameter values from the parameter file
    // if (useprm) readprm();
}
