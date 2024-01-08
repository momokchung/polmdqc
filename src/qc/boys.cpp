// Author: Moses KJ Chung
// Year:   2023

#include "boys.h"
#include "init.h"
#include "mathUtils.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace boys
{
///////////////////////////////
//                           //
//  boys  --  Boys function  //
//                           //
///////////////////////////////

constexpr int M = 52;
constexpr int intervalLength = 67999;
// constexpr int N = intervalLength * 4 * (M + 1);

char* memblock;
real* coefficients;

//////////////////////////////////////////////
//                                          //
//  initBoys  --  initialize Boys function  //
//                                          //
//////////////////////////////////////////////

void initBoys()
{
    std::streampos size;
    std::string fileName = init::cwd;
    fileName.append("/f0052.bin");
    std::ifstream file (fileName, std::ios::in|std::ios::binary|std::ios::ate);
    size = file.tellg();
    // std::cout << "size = " << size << "\n";

    memblock = new char [size];
    file.seekg (0, std::ios::beg);
    file.read (memblock, size);
    file.close();
    coefficients = (real*)memblock;
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//  boysIntegralPoly  --  returns Chebyshev polynomial approximation to Boys function  //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

void boysIntegralPoly(real t, int m, std::vector<real>& boysPoly)
{
    constexpr real d = 5e-4;
    real x = t/d;
    if (t >= 0. and t < 34)
    {
        int j = int(x);
        for (int i = 0; i <= m; ++i)
        {
            int index = intervalLength * 4 * i + j * 4;
            boysPoly[i] = coefficients[index] + x * (coefficients[index + 1] + x * (coefficients[index + 2] + x * coefficients[index + 3]));
        }
    }
    else
    {
        real t2 = t * 2.;
        boysPoly[0] = sqrt(unitsqm::pi / 2. / t2);
        for (int i = 1; i <= m; ++i)
        {
            // boysPoly[i] = boysPoly[i] * (2 * i - 1) / t2; // less accurate, fast
            boysPoly[i] = (boysPoly[i-1] * (2 * i - 1) - exp(-t))/ t2; // more accurate, slow
        }
    }
    return ;
}

//////////////////////////////////////////////
//                                          //
//  cleanupBoys  --  cleanup Boys function  //
//                                          //
//////////////////////////////////////////////

void cleanupBoys()
{
    delete[] memblock;
    return;
}
}
