// Author: Moses KJ Chung
// Year:   2024

#include "alfp.h"
#include "gettext.h"
#include "kalf.h"
#include "keys.h"
#include "openmp.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  kalf  --  AlphaMol parameter assignment  //
//                                           //
///////////////////////////////////////////////

// "kalf" assigns the parameters to be used in computing
// surface area and volume using AlphaMol/AlphaMol2

void kalf()
{
    int next;
    std::string str;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // set default parameters
    alfmeth = AlfMethod::AlphaMol2;
    alfsort = AlfSort::KDTree;
    nthdDlny = nthread;

    // process keywords containing AlphaMol parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "ALF-SORT") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            if (!(iss >> str)) continue;
            if (str == "NONE") {
                alfsort = AlfSort::None;
            }
            else if (str == "SORT3D") {
                alfsort = AlfSort::Sort3D;
            }
            else if (str == "BRIO") {
                alfsort = AlfSort::BRIO;
            }
            else if (str == "SPLIT") {
                alfsort = AlfSort::Split;
            }
            else if (str == "KDTREE") {
                alfsort = AlfSort::KDTree;
            }
        }
        else if (keyword == "ALF-METHOD") {
            string = record.substr(next);
            upcase(string);
            iss.clear();
            iss.str(string);
            if (!(iss >> str)) continue;
            if (str == "ALPHAMOL") {
                alfmeth = AlfMethod::AlphaMol;
            }
            else if (str == "ALPHAMOL2") {
                alfmeth = AlfMethod::AlphaMol2;
            }
        }
        if (keyword == "DLNY-THREADS") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> nthdDlny)) continue;
            if (nthdDlny < 1) nthdDlny = 1;
        }
    }
}
}
