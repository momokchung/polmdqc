// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////
//                                          //
//  alfc  --  Alpha shape theory variables  //
//                                          //
//////////////////////////////////////////////

// alfeps   epsilon value

constexpr real alfeps = 1e-5;

constexpr int other3[4][3] = {
    {1, 2, 3},
    {0, 2, 3},
    {0, 1, 3},
    {0, 1, 2}
};

constexpr int face_info[6][2] = {
    {0, 1},
    {0, 2},
    {0, 3},
    {1, 2},
    {1, 3},
    {2, 3}
};

constexpr int face_pos[6][2] = {
    {1, 0},
    {2, 0},
    {3, 0},
    {2, 1},
    {3, 1},
    {3, 2}
};

constexpr int pair[6][2] = {
    {2, 3},
    {1, 3},
    {1, 2},
    {0, 3},
    {0, 2},
    {0, 1}
};
}
