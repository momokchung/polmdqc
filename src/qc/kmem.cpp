// Author: Moses KJ Chung
// Year:   2024

#include "gettext.h"
#include "keys.h"
#include "kmem.h"
#include "mem.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////
//                                       //
//  kmem  --  assign memory information  //
//                                       //
///////////////////////////////////////////

// "kmem" assigns memory information

void kmem()
{
    int next;
    int mem;
    std::string memunit;
    std::string keyword,record,string;
    std::istringstream iss;

    // process keywords containing memory parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "MEMORY") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            mem = 0;
            memunit = "";
            iss >> mem >> memunit;
            if (mem <= 0) continue;
            if (memunit == "MB") memory = mem;
            else if (memunit == "GB") memory = mem * 1024;
        }
    }
}
}
