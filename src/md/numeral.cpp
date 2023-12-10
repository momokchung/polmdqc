// Author: Moses KJ Chung
// Year:   2023

#include "numeral.h"
#include <iomanip>
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  numeral  --  convert number to text string  //
//                                              //
//////////////////////////////////////////////////

// "numeral" converts an input integer number into the
// corresponding text numeral padded with zeros
// 
// number  integer value of the number to be transformed
// width    total length of string

std::string numeral(int number, int width)
{
    std::string numberString = std::to_string(number);
    if (numberString.length() > width)
        return numberString;

    // Create an output string stream
    std::ostringstream oss;

    // Use std::setw and std::setfill to pad the number with zeros
    oss << std::setw(width) << std::setfill('0') << number;

    // Convert the output stream to a string
    std::string paddedNumber = oss.str();

    return paddedNumber;
}
}
