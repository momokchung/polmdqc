// Author: Moses KJ Chung
// Year:   2024

#include "alfp.h"
#include "gettext.h"
#include "kalf.h"
#include "keys.h"
#include "openmp.h"
#include "upcase.h"
#include <cmath>
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

inline int roundDownToPowerOf2(int x) {
    if (x <= 0)
        return 1;
    int exponent = static_cast<int>(std::floor(std::log2(x)));
    return static_cast<int>(std::pow(2, exponent));
}

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
    alfdigit = 8;
    alfnthd = roundDownToPowerOf2(nthread);
    alfsos = true;

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
        else if (keyword == "DLNY-THREADS") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> alfnthd)) continue;
            if (alfnthd < 1) alfnthd = 1;
        }
        else if (keyword == "ALF-DIGITS") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> alfdigit)) continue;
            if (alfdigit < 8) alfdigit = 8;
            // round alfdigit down to nearest even integer
            if ((alfdigit % 2) != 0) alfdigit -= 1;
        }
        else if (keyword == "ALF-NOSOS") {
            alfsos = false;
        }
    }
}
}
