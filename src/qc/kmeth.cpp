// Author: Moses KJ Chung
// Year:   2024

#include "gettext.h"
#include "keys.h"
#include "kmeth.h"
#include "methqm.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  kmeth  --  assign QM theory/method parameters  //
//                                                 //
/////////////////////////////////////////////////////

// "kmeth" assigns QM theory/method parameters

void kmeth()
{
    int next;
    std::string meth;
    std::string keyword,record,string;
    std::istringstream iss;

    // process keywords containing theory/method parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "METHOD") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            meth = "";
            iss >> meth;
            if (meth == "HF") qmmethod = QMMethodType::HF;
            else if (meth == "MP2") qmmethod = QMMethodType::MP2;
            else if (meth == "CCSD") qmmethod = QMMethodType::CCSD;
            else if (meth == "CCSD(T)") qmmethod = QMMethodType::CCSDt;
            else if (meth == "CCSDT") qmmethod = QMMethodType::CCSDT;
            else if (meth == "SAPT0") qmmethod = QMMethodType::SAPT0;
            else if (meth == "SAPT2") qmmethod = QMMethodType::SAPT2;
            else if (meth == "SAPT2+") qmmethod = QMMethodType::SAPT2p;
            else if (meth == "ALMO") qmmethod = QMMethodType::ALMO;
        }
    }
}
}
