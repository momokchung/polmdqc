///////////////////////////////////////////////////////
//                                                   //
//  torphase.cpp  --  torsional amplitude and phase  //
//                                                   //
///////////////////////////////////////////////////////

// "torphase" sets the n-fold amplitude and phase values
// for each torsion via sorting of the input parameters


#include "torphase.h"

void torphase(int (&ft)[6], double (&vt)[6], double (&st)[6])
{
    double ampli[6];
    double phase[6];

    // copy the input fold, amplitude and phase angles
    for (int i = 0; i < 6; i++){
        ampli[i] = vt[i];
        phase[i] = st[i];
        vt[i] = 0.;
        st[i] = 0.;
    }

    // shift the phase angles into the standard range
    for (int i = 0; i < 6; i++) {
        while (phase[i] < -180.) {
            phase[i] += 360.;
        }
        while (phase[i] > 180.) {
            phase[i] -= 360.;
        }
    }

    // convert input torsional parameters to storage format
    for (int i = 0; i < 6; i++) {
        int k = ft[i];
        if (k>=0 and k<=5) {
            vt[k] = ampli[i];
            st[k] = phase[i];
        }
    }
}
