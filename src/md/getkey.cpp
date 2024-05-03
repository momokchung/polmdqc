// Author: Moses KJ Chung
// Year:   2023

#include "argue.h"
#include "fatal.h"
#include "files.h"
#include "getnumb.h"
#include "gettext.h"
#include "inquire.h"
#include "keys.h"
#include "openmp.h"
#include "suffix.h"
#include "upcase.h"
#include "version.h"
#include <fstream>
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  getkey  --  find and store contents of keyfile  //
//                                                  //
//////////////////////////////////////////////////////

// "getkey" finds a valid keyfile and stores its contents as
// line images for subsequent keyword parameter searching

void getkey()
{
    int next;
    bool exist,header;
    std::string keyword;
    std::string keyfile;
    std::string record;
    std::string string;
    std::istringstream iss;

    // check for a keyfile specified on command line
    exist = false;
    for (int i = 0; i < narg - 1; i++) {
        string = arg[i];
        upcase(string);
        if (string.substr(0, 2) == "-K") {
            keyfile = arg[i + 1];
            suffix (keyfile, "key", "old");
            exist = inquireFile(keyfile);
            if (!exist) {
                printf("\n GETKEY  --  Keyfile Specified on Command Line was not Found\n");
                fatal();
            }
        }
    }

    // try to get keyfile from base name of current system
    if (!exist) {
        keyfile = filename.substr(0,leng) + ".key";
        version(keyfile, "old");
        exist = inquireFile(keyfile);
    }

    // check for the existence of a generic keyfile
    if (!exist) {
        if (ldir == 0) {
            keyfile = "polmdqc.key";
        }
        else {
            keyfile = filename.substr(0,ldir) + "polmdqc.key";
        }
        version(keyfile, "old");
        exist = inquireFile(keyfile);
    }

    // read the keyfile and store it for latter use
    nkey = 0;
    if (exist) {
        std::ifstream file(keyfile);
        while (std::getline(file, record)) {
            keyline[nkey] = record;
            nkey++;
            if (nkey >= maxkey) {
                printf("\n GETKEY  --  Keyfile Too Large; Increase MAXKEY\n");
                fatal();
            }
        }
        file.close();
    }

    // convert underbar characters to dashes in all keywords
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record, keyword, next);
        for (int j = 0; j < next; j++)
        {
            if (record[j] == '_') record[j] = '-';
        }
        keyline[i] = record;
    }

    // check for comment lines to be echoed to the output
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        if (keyword == "ECHO") {
            if (header) {
               header = false;
               printf("\n");
            }
            printf("%s\n", string.c_str());
        }
    }

    // set number of OpenMP threads for parallelization
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        upcase(record);
        gettext(record,keyword,next);
        string = record.substr(next);
        if (keyword == "OPENMP-THREADS") {
            iss.clear();
            iss.str(string);
            if (!(iss >> nthread)) break;
            if (nthread == 0) nthread = 1;
            omp_set_num_threads(nthread);
        }
    }

    // check for number of OpenMP threads on command line
    for (int i = 0; i < narg - 1; i++) {
        string = arg[i];
        upcase (string);
        if (string.substr(0, 2) == "-T") {
            next = 0;
            string = arg[i + 1];
            getnumb(string,nthread,next);
            if (nthread == 0) nthread = 1;
            omp_set_num_threads(nthread);
        }
    }
}
}
