// Author: Moses KJ Chung
// Year:   2024

#include "gettext.h"
#include "keys.h"
#include "ksym.h"
#include "sym.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////
//                                        //
//  ksym  --  assign symmetry parameters  //
//                                        //
////////////////////////////////////////////

// "ksym" assigns symmetry to current structure

void ksym()
{
    int next;
    std::string sym;
    std::string keyword,record,string;
    std::istringstream iss;

    // implement method to detect symmetry

    // process keywords containing symmetry parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "SYMMETRY") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            sym = "";
            iss >> sym;
            if (sym == "C1") symmetry = Symmetry::C1;
            else if (sym == "CI") symmetry = Symmetry::Ci;
            else if (sym == "C2") symmetry = Symmetry::C2;
            else if (sym == "CS") symmetry = Symmetry::Cs;
            else if (sym == "D2") symmetry = Symmetry::D2;
            else if (sym == "C2V") symmetry = Symmetry::C2v;
            else if (sym == "C2H") symmetry = Symmetry::C2h;
            else if (sym == "D2H") symmetry = Symmetry::D2h;
        }
    }
    if (symmetry == Symmetry::C1) printf("C1\n");
    else if (symmetry == Symmetry::Ci) printf("Ci\n");
    else if (symmetry == Symmetry::C2) printf("C2\n");
    else if (symmetry == Symmetry::Cs) printf("Cs\n");
    else if (symmetry == Symmetry::D2) printf("D2\n");
    else if (symmetry == Symmetry::C2v) printf("C2v\n");
    else if (symmetry == Symmetry::C2h) printf("C2h\n");
    else if (symmetry == Symmetry::D2h) printf("D2h\n");
}
}
