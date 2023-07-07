/////////////////////////////////////////////////////////
//                                                     //
//  basefile.cpp  --  get base prefix from a filename  //
//                                                     //
/////////////////////////////////////////////////////////

// "basefile" extracts from an input filename the portion
// consisting of any directory name and the base filename;
// also reads any keyfile and sets information level values


#include "basefile.h"
#include "control.h"
#include "files.h"
#include "getkey.h"

void basefile(std::string& string)
{
    // account for home directory abbreviation in filename
    if (string.substr(0,2).compare("~/") == 0) {
        char* homeDir = getenv("HOME");
        if (homeDir != nullptr) {
            std::string prefix = std::string(homeDir);
            string = prefix + string.substr(1).c_str();
        }
    }
    
    // store the input filename and find its full length
    filename = string;
    leng = string.length();
    
    // count the number of characters prior to any extension
    int k = leng;
    for (int i = 0; i < leng; i++) {
        char letter = string[i];
        if (letter == '/') k = leng;
        else if (letter == '\\') k = leng;
        else if (letter == ']') k = leng;
        else if (letter == ':') k = leng;
        else if (letter == '~') k = leng;
        else if (letter == '.') k = i;
    }

    // set the length of the base file name without extension
    leng = std::min(leng,k);

    // find the length of any directory name prefix
    k = 0;
    for (int i = leng; i > 0; i--) {
        char letter = string[i-1];
        if (letter == '/') k = i;
        else if (letter == '\\') k = i;
        else if (letter == ']') k = i;
        else if (letter == ':') k = i;
        else if (letter == '~') k = i;
        else if (letter == '.') k = i;
        if (k != 0) break;
    }
    ldir = k;

    // read and store the keywords from the keyfile
    getkey();

    // get the information level and output style
    control();
}
