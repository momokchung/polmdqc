// Author: Moses KJ Chung
// Year:   2023

#include "ptable.h"
#include <iostream>
#include <string>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  initatom  --  setup atoms in periodic table  //
//                                               //
///////////////////////////////////////////////////

// "initatom" sets the atomic symbol, standard atomic weight,
// van der Waals radius and covalent radius for each element in
// the periodic table
// 
// literature references:
// 
// J. Emsley, The Elements, 3rd Edition, Oxford University Press,
// (1999)  [relative atomic masses]
// 
// J. Meija, T. B. Coplen, M. Berglund, W. A. Brand, P. De Bievre,
// M. Groning, N. E. Holden, J. Irrgeher, R. D. Loss, T. Walczyk and
// R. Prohaska, Atomic Weights of the Elements 2013, Pure and Applied
// Chemistry, 88, 265-291 (2016)  [standard atomic weights]
// 
// A. Bondi, van der Waals Volumes and Radii, Journal of Physical
// Chemistry, 68, 441-451 (1964)  [original vdw radii; not used]
// 
// S. Alvarez, "A Cartography of the van der Waals Territories",
// Dalton Transactions, 42, 8617-8636 (2013)  [vdw radii for most
// elements 1-99]
// 
// B. Cordero, V. Gomez. A. E. Platero-Prats, M. Reves,
// J. Echeverria, E. Cremades, F. Barragan and S. Alverez,
// "Covalent Radii Revisited", Dalton Transactions, 2832-2838 (2008)
// [covalent radii for elements 1-96]
// 
// P. Pyykko and M. Atsumi, "Molecular Single-Bond Covalent Radii
// for Elements 1-118", Chemistry- A European Journal, 15, 187-197
// (2009)  [covalent radii for elements 97-112]

void initatom()
{
    // atomic symbol for each element
    std::string asym[maxele] = {
        "H  ", "He ", "Li ", "Be ", "B  ", "C  ", "N  ",
        "O  ", "F  ", "Ne ", "Na ", "Mg ", "Al ", "Si ",
        "P  ", "S  ", "Cl ", "Ar ", "K  ", "Ca ", "Sc ",
        "Ti ", "V  ", "Cr ", "Mn ", "Fe ", "Co ", "Ni ",
        "Cu ", "Zn ", "Ga ", "Ge ", "As ", "Se ", "Br ",
        "Kr ", "Rb ", "Sr ", "Y  ", "Zr ", "Nb ", "Mo ",
        "Tc ", "Ru ", "Rh ", "Pd ", "Ag ", "Cd ", "In ",
        "Sn ", "Sb ", "Te ", "I  ", "Xe ", "Cs ", "Ba ",
        "La ", "Ce ", "Pr ", "Nd ", "Pm ", "Sm ", "Eu ",
        "Gd ", "Tb ", "Dy ", "Ho ", "Er ", "Tm ", "Yb ",
        "Lu ", "Hf ", "Ta ", "W  ", "Re ", "Os ", "Ir ",
        "Pt ", "Au ", "Hg ", "Tl ", "Pb ", "Bi ", "Po ",
        "At ", "Rn ", "Fr ", "Ra ", "Ac ", "Th ", "Pa ",
        "U  ", "Np ", "Pu ", "Am ", "Cm ", "Bk ", "Cf ",
        "Es ", "Fm ", "Md ", "No ", "Lr ", "Rf ", "Db ",
        "Sg ", "Bh ", "Hs ", "Mt ", "Ds ", "Rg ", "Cn "
    };

    // standard atomic weight for each element
    double amas[maxele] = {
        1.008,   4.003,   6.941,   9.012,  10.811,
        12.011,  14.007,  15.999,  18.998,  20.180,
        22.990,  24.305,  26.982,  28.086,  30.974,
        32.066,  35.453,  39.948,  39.098,  40.078,
        44.956,  47.867,  50.942,  51.996,  54.938,
        55.845,  58.933,  58.693,  63.546,  65.380,
        69.723,  72.630,  74.922,  78.971,  79.904,
        83.798,  85.468,  87.620,  88.906,  91.224,
        92.906,  95.950,  98.906, 101.070, 102.910,
        106.420, 107.870, 112.410, 114.820, 118.710,
        121.760, 127.600, 126.900, 131.290, 132.910,
        137.330, 138.910, 140.120, 140.910, 144.240,
        144.913, 150.360, 151.960, 157.250, 158.930,
        162.500, 164.930, 167.260, 168.930, 173.050,
        174.970, 178.490, 180.950, 183.840, 186.210,
        190.230, 192.220, 195.080, 196.970, 200.590,
        204.383, 207.200, 208.980, 208.982, 209.987,
        222.017, 223.020, 226.025, 227.027, 232.038,
        231.036, 238.029, 237.048, 244.064, 243.061,
        247.070, 247.070, 251.080, 252.083, 257.095,
        258.098, 259.101, 262.110, 267.122, 270.131,
        269.129, 270.133, 270.134, 278.156, 281.165,
        281.166, 285.177
    };

    // van der Waals radius for each element (Angstroms)
    double vrad[maxele] = {
        1.20, 1.43, 2.12, 1.98, 1.91, 1.77,
        1.66, 1.50, 1.46, 1.58, 2.50, 2.51,
        2.25, 2.19, 1.90, 1.89, 1.82, 1.83,
        2.73, 2.62, 2.58, 2.46, 2.42, 2.45,
        2.45, 2.44, 2.40, 2.40, 2.38, 2.39,
        2.32, 2.29, 1.88, 1.82, 1.86, 2.25,
        3.21, 2.84, 2.75, 2.52, 2.56, 2.45,
        2.44, 2.46, 2.44, 2.15, 2.53, 2.49,
        2.43, 2.42, 2.47, 1.99, 2.04, 2.06,
        3.48, 3.03, 2.98, 2.88, 2.92, 2.95,
        2.90, 2.90, 2.87, 2.83, 2.79, 2.87,
        2.81, 2.83, 2.79, 2.80, 2.74, 2.63,
        2.53, 2.57, 2.49, 2.48, 2.41, 2.29,
        2.32, 2.45, 2.47, 2.60, 2.54, 2.93,
        2.88, 2.71, 2.82, 2.81, 2.80, 2.93,
        2.88, 2.71, 2.82, 2.81, 2.83, 3.05,
        3.40, 3.05, 2.70, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.00
    };

    // covalent radius for each element (Angstroms)
    double crad[maxele] = {
        0.31, 0.28, 1.28, 0.96, 0.84, 0.76,
        0.71, 0.66, 0.57, 0.58, 1.66, 1.41,
        1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
        2.03, 1.76, 1.70, 1.60, 1.53, 1.39,
        1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
        1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
        2.20, 1.95, 1.90, 1.75, 1.64, 1.54,
        1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
        1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
        2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
        1.99, 1.98, 1.98, 1.96, 1.94, 1.92,
        1.92, 1.89, 1.90, 1.87, 1.87, 1.75,
        1.70, 1.62, 1.51, 1.44, 1.41, 1.36,
        1.36, 1.32, 1.45, 1.46, 1.48, 1.40,
        1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
        2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
        1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
        1.61, 1.57, 1.49, 1.43, 1.41, 1.34,
        1.29, 1.28, 1.21, 1.22
    };

    // set the symbol, weight and radii for each element
    for (int i = 0; i < maxele; i++) {
        atmass[i] = amas[i];
        elemnt[i] = asym[i];
        if (vrad[i] == 0.0) vrad[i] = 2.0;
        vdwrad[i] = vrad[i];
        covrad[i] = crad[i];
    }
}
}
