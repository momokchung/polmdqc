///////////////////////////////////////////////////////
//                                                   //
//  version.cpp  --  create version number for file  //
//                                                   //
///////////////////////////////////////////////////////

// "version" checks the name of a file about to be opened; if
// if "old" status is passed, the name of the highest current
// version is returned; if "new" status is passed the filename
// of the next available unused version is generated


#include "inquire.h"
#include "lowcase.h"
#include "nextarg.h"
#include "output.h"
#include "version.h"
#include <iostream>

void version(std::string& string, std::string status)
{
    // process the filename and status variables
    lowcase(status);
    int leng = string.length();

    // no change is needed if the file doesn't exist
    bool exist = false;
    if (leng != 0)  exist = inquire(string);
    if (!exist)  return;

    // set initial values for the current and next versions
    std::string newfile = string;
    std::string oldfile = string;

    // append an artificial version number to the filename;
    // currently handles up to 10000 versions of a file
    if (!noversion) {
        int i = 1;
        while (exist) {
            i++;
            oldfile = newfile;
            int thousand = i / 1000;
            int hundred = (i - 1000 * thousand) / 100;
            int tens = (i - 1000 * thousand - 100 * hundred) / 10;
            int ones = i - 1000 * thousand - 100 * hundred - 10 * tens;
            std::string number;
            if (thousand != 0) {
                number = std::to_string(thousand)
                       + std::to_string(hundred)
                       + std::to_string(tens)
                       + std::to_string(ones);
            } else if (hundred != 0) {
                number = std::to_string(hundred)
                       + std::to_string(tens)
                       + std::to_string(ones);
            } else if (tens != 0) {
                number = std::to_string(tens)
                       + std::to_string(ones);
            } else {
                number = std::to_string(ones);
            }
            newfile = string + "_" + number;
            exist = inquire(newfile);
        }
    }

    // set the file name based on the requested status
    if (status == "old") {
        string = oldfile;
    } else if (status == "new") {
        string = newfile;
        exist = inquire(string);
        if (exist) {
            nextarg(string, exist);
            if (exist) {
                exist = inquire(string);
            } else {
                exist = true;
            }
            while (exist) {
                printf("\n Enter File Name for Coordinate Output:  ");
                std::getline(std::cin, string);
                exist = inquire(string);
            }
        }
    }
}
