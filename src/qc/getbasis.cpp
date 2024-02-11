// Author: Moses KJ Chung
// Year:   2024

#include "basis.h"
#include "fatal.h"
#include "getbasis.h"
#include "inquire.h"
#include "keys.h"
#include "lowcase.h"
#include "readbasis.h"
#include "rootdir.h"
#include "setbasis.h"
#include "suffix.h"
#include "upcase.h"
#include <fstream>
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  getbasis  --  get basis set parameter file  //
//                                              //
//////////////////////////////////////////////////

// "getbasis" finds the basis set parameter file
// and then opens and reads the basis functions

void getbasis()
{
    int next;
    bool exist;
    std::string gbsfile;
    std::string gbsname;
    std::string keyword;
    std::string record;
    std::istringstream iss;

    // search the keyword list for the basis set name
    for (int i = 0; i < nkey; i++) {
        next = 1;
        record = keyline[i];
        iss.clear();
        keyword = "";
        iss.str(record);
        iss >> keyword;
        upcase(keyword);
        if (keyword == "BASISSET" or keyword == "BASIS") {
            iss >> gbsname;
        }
    }

    // set default basis set name
    upcase(gbsname);
    bssname = gbsname;

    // add .gbs extension
    lowcase(gbsname);
    gbsname.append(".gbs");

    // append gbsname to gbsfile and see if file exists
    gbsfile = rootdir();
    gbsfile.append("/basis/");
    gbsfile.append(gbsname);
    exist = inquireFile(gbsfile);

    if (!exist) {
        printf("\n GETBASIS  --  No Basis Set Given; Specify Basis Set\n");
        fatal();
    }

    // read the basis set file and store it for latter use
    nbss = 0;
    std::ifstream file(gbsfile);
    while (std::getline(file, record)) {
        bssline[nbss] = record;
        nbss++;
        if (nbss >= maxbss) {
            printf("\n GETBASIS  --  Basis Set File Too Large; Increase MAXBSS\n");
            fatal();
        }
    }
    file.close();

    // allocate memory for basis set values
    setbasis();

    // get basis set values from the basis set file
    readbasis();
}
}
