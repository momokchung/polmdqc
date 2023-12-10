// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  ksolv  --  solvation parameter assignment  //
//                                             //
/////////////////////////////////////////////////

void ksolv();

void kgb();
void kgk();
void khpmf();
void knp();
void kpb();
void ksa();
void setrad(std::string radtyp);
}
