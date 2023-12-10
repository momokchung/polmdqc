// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <Eigen/Dense>
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  rotpole  --  rotate multipoles to global frame  //
//                                                  //
//////////////////////////////////////////////////////

void rotpole(std::string rotMode);

void rotrpole(std::string rotMode);

void rotmat(int i, Eigen::Matrix<double, 3, 3>& a, bool& planar);

void rotsite(int ii, Eigen::Matrix<double, 3, 3>& a, bool& planar, std::vector<std::vector<double>>& pole, std::vector<std::vector<double>>& rpole);
}
