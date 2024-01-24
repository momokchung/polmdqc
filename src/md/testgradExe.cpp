// Author: Moses KJ Chung
// Year:   2024

#include "inform.h"
#include "testgrad.h"

///////////////////
//               //
//  testgradExe  //
//               //
///////////////////

int main(int argc, char** argv)
{
    polmdqc::test = false;
    polmdqc::testgrad(argc, argv);
}
