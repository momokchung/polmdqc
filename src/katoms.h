///////////////////////////////////
//                               //
//  katoms.h  --  atoms in file  //
//                               //
///////////////////////////////////


#pragma once
#include "init.h"
#include <map>
#include <string>
#include <vector>

namespace atoms 
{
enum class LengthUnit
{
    angstrom,
    bohr
};

enum class Symmetry
{
    C1,
    C2,
    Ci,
    Cs,
    D2,
    C2h,
    C2v,
    D2h
};

extern int n;
extern int nElec;
extern int charge;
extern int multiplicity;
extern bool com;
extern bool reorient;
extern std::vector<int> core;
extern std::vector<real> coordx;
extern std::vector<real> coordy;
extern std::vector<real> coordz;
extern std::string name;
extern std::vector<std::string> atom;
extern LengthUnit lengthUnit;
extern Symmetry symmetry;

extern std::map<std::string, int> atomicNumber;
extern std::map<std::string, real> atomicWeight;

void readxyz(std::string fileName);
}
