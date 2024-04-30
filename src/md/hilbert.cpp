// Author: Moses KJ Chung
// Year:   2024

#include "hilbert.h"
#include "hilbrt.h"
#include <stdio.h>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  initHilbert  --  initialize the permutation table  //
//                                                     //
/////////////////////////////////////////////////////////

// "initHilbert" initializes the permutation table.
//
// The table 'transgc' has 8 x 3 x 8 entries. It contains
// all possible Gray code sequences traveled by the 1st
// order Hilbert curve in 3 dimensions. The first column
// is the Gray code of the entry point of the curve, and
// the second column is the direction (0, 1, or 2, 0 means
// the x-axis) where the exit point of curve lies.
// The table 'tsb1mod3' contains the numbers of trailing
// set '1' bits of the indices from 0 to 7, modulo by '3'.
// The code for generating this table is from:
// https://graphics.stanford.edu/~seander/bithacks.html

void initHilbert(int ndim)
{
    int gc[8],N,mask,travel_bit;
    int c,f,g,k,v;

    for (int i = 0; i < 8; i++) gc[i] = 0;

    N = (ndim == 2) ? 4 : 8;
    mask = (ndim == 2) ? 3 : 7;

    // Generate the Gray code sequence.
    for (int i = 0; i < N; i++) {
        gc[i] = i ^ (i >> 1);
    }

    for (int e = 0; e < N; e++) {
        for (int d = 0; d < ndim; d++) {
            // Calculate the end point (f).
            f = e ^ (1 << d); // Toggle the d-th bit of 'e'.
            // travel_bit = 2**p, the bit we want to travel. 
            travel_bit = e ^ f;
            for (int i = 0; i < N; i++) {
                // // Rotate gc[i] left by (p + 1) % n bits.
                k = gc[i] * (travel_bit * 2);
                g = ((k | (k / N)) & mask);
                // Calculate the permuted Gray code by xor with the start point (e).
                transgc[e][d][i] = (g ^ e);
            }
        } // d
    } // e

    // Count the consecutive '1' bits (trailing) on the right.
    tsb1mod3[0] = 0;
    for (int i = 1; i < N; i++) {
        v = ~i; // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        for (c = 0; v; c++) {
            v >>= 1;
        }
        tsb1mod3[i] = c % ndim;
    }
}
}
