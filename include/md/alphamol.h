/* ===============================================================================================
   AlphaMol: a program for computing geometric measures of a union of balls

   Author:  Patrice Koehl  (collaboration with Herbert Edelsbrunner
   Date:    9/22/2019
   Version: 1
   =============================================================================================== */

#pragma once
#include "alfcx.h"
#include "delcx.h"
#include "edge.h"
#include "face.h"
#include "tetrahedron.h"
#include "vertex.h"
#include "volumes.h"
#include <bitset>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

void alphamol(double r_h2o, int flag_deriv);
}
