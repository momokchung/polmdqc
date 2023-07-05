///////////////////////////////////////////////////////
//                                                   //
//  suffix.cpp  --  test for default file extension  //
//                                                   //
///////////////////////////////////////////////////////

// "suffix" checks a filename for the presence of an extension,
// and appends an extension and version if none is found


#include "inquire.h"
#include "suffix.h"
#include "version.h"

void suffix(std::string& string, std::string extension, std::string status)
{
    // get the full length of the current filename
    int leng = string.length();
    int lext = extension.length();

    // check for an extension on the current filename
    int k = leng;
    for (int i = 0; i < leng; i++) {
        char letter = string[i];
        if (letter == '/')
            k = leng;
        if (letter == '\\')
            k = leng;
        if (letter == ']')
            k = leng;
        if (letter == ':')
            k = leng;
        if (letter == '~')
            k = leng;
        if (letter == '.')
            k = i;
    }

    // append an extension or version as appropriate
    if (k == leng) {
        bool exist = false;
        if (leng != 0) {
            exist = inquire(string);
        }
        if (!exist) {
            string = string + "." + extension;
            version(string, status);
        }
    }
    else if (status == "new") {
        version(string, status);
    }
}
