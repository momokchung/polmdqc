// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <Eigen/Dense>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  rotpole  --  rotate multipoles to global frame  //
//                                                  //
//////////////////////////////////////////////////////

enum class RotMode
{
    None,
    Mpole,
    Repel,
};

void rotpole(RotMode rotMode);

void rotrpole(RotMode rotMode);

void rotmat(int i, Eigen::Matrix<real, 3, 3>& a, bool& planar);

void rotsite(int ii, Eigen::Matrix<real, 3, 3>& a, bool& planar, std::vector<std::vector<real>>& pole, std::vector<std::vector<real>>& rpole);
}
