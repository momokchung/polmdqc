// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////
//                                              //
//  prmkey  --  interpret force field keywords  //
//                                              //
//////////////////////////////////////////////////

void prmkey(std::string record);

void potoff();

void valoff();

void nbondoff();
}
