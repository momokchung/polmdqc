/////////////////////////////////
//                             //
//  boys.h  --  Boys function  //
//                             //
/////////////////////////////////


#pragma once
#include "init.h"
#include "units.h"
#include <vector>

namespace boys
{
extern char* memblock;
extern real* coefficients;

void initBoys();
void boysIntegralPoly(real t, int m, std::vector<real>& boysPoly);
void cleanupBoys();
}
