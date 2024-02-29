// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "cheb.h"
#include <cmath>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  boysref  --  Boys function reference evaluation  //
//                                                   //
///////////////////////////////////////////////////////

// "boysref" evaluates the reference Boys function
//
// literature reference:
//
// I. Shavitt, "The Gaussian Function in Calculations
// of Statistical Mechanics and Quantum Mechanics",
// in: Methods in Computational Physics, 2, 1-45, (1963).
//
// E. F. Valeev, "Libint: A library for the evaluationof 
// molecular integrals of many-body operators over Gaussian
// functions", Version 2.7.2, http://libint.valeyev.net/

inline void boysref(double* Fm, int mmax, double t)
{
    if (t < chebtmax) {
        double et = std::exp(-t);
        for (int m = 0; m <= mmax; m++) {
            constexpr double half = double(1)/2;
            double denom = (m + half);
            double term = et / (2 * denom);
            double old_term = 0;
            double sum = term;
            const double epsilon = std::numeric_limits<double>::epsilon();
            const double epsilon_divided_10 = epsilon / 10;
            do {
                denom += 1;
                old_term = term;
                term = old_term * t / denom;
                sum += term;
                // rel_error = term / sum , hence iterate until rel_error = epsilon
                // however, must ensure that contributions are decreasing to ensure that omitted contributions are smaller than epsilon
            } while (term > sum * epsilon_divided_10 || old_term < term);

            Fm[m] = sum;
        }
    }
    else {
        constexpr double two_over_sqrt_pi{1.128379167095512573896158903121545171688101258657997713688171443421284936882986828973487320404214727};
        constexpr double K = 1/two_over_sqrt_pi;

        double t2 = 2*t;
        double et = std::exp(-t);
        double sqrt_t = sqrt(t);
        Fm[0] = K * std::erf(sqrt_t) / sqrt_t;
        for (int m = 0; m <= mmax-1; m++) {
            Fm[m+1] = ((2*m + 1)*Fm[m] - et)/(t2);
        }
    }
}
}
