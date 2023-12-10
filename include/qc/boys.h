// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include "unitsqm.h"
#include <vector>

namespace boys
{
///////////////////////////////
//                           //
//  boys  --  Boys function  //
//                           //
///////////////////////////////

extern char* memblock;
extern real* coefficients;

void initBoys();
void boysIntegralPoly(real t, int m, std::vector<real>& boysPoly);
void cleanupBoys();
}
