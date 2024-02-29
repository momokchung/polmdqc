// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "cheb.h"
#include <cmath>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  boyscheb  --  Boys function Chebyshev evaluation  //
//                                                    //
////////////////////////////////////////////////////////

// "boyscheb" evaluates the Boys function using the 
// Chebyshev expansion
//
// literature reference:
//
// P. M. W. Gill, B. G. Johnson, and J. A. Pople, "Two-Electron
// Repulsion Integrals Over Gaussian s Functions", International
// Journal of Quantum Chemistry, 40, 745-752, (1991).
//
// E. F. Valeev, "Libint: A library for the evaluationof 
// molecular integrals of many-body operators over Gaussian
// functions", Version 2.7.2, http://libint.valeyev.net/

inline void boyscheb(double* Fm, int mmax, double x)
{
    if (x > chebtmax) {
        const double one_over_x = 1/x;

        // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
        Fm[0] = 0.88622692545275801365 * std::sqrt(one_over_x);
        if (mmax == 0)
            return;

        // // see Eq. (9.8.13) in Helgaker-Jorgensen-Olsen
        // // FAST version; omits "-e^(-x)/(2x)" term
        // for (int i = 1; i <= mmax; i++) {
        //     Fm[i] = Fm[i - 1] * ihalf[i] * one_over_x;
        // }

        // SLOW version; includes "-e^(-x)/(2x)" term
        double ex2 = std::exp(-x)/2;
        for (int i = 1; i <= mmax; i++) {
            Fm[i] = (Fm[i - 1] * ihalf[i] - ex2) * one_over_x;
        }
        return;
    }

    // set up Chebyshev interpolation
    const int mmin = 0;
    int iv = int(x/chebdlta);
    double t = x/chebdlta2 - 2*iv - 1;

    // compute Chebyshev interpolation
    for (int m = mmin; m <= mmax; m++) {
        int index = m*(intorder+1);
        Fm[m] = chebtble[iv][index + 0]
            + t * (chebtble[iv][index + 1]
            + t * (chebtble[iv][index + 2]
            + t * (chebtble[iv][index + 3]
            + t * (chebtble[iv][index + 4]
            + t * (chebtble[iv][index + 5]
            + t * (chebtble[iv][index + 6]
            + t * (chebtble[iv][index + 7])))))));
    }
}
}
