// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  fft  --  Fast Fourier transform control values  //
//                                                  //
//////////////////////////////////////////////////////

// maxprime   maximum number of prime factors of FFT dimension
// iprime     prime factorization of each FFT dimension (FFTPACK)
// planf      pointer to forward transform data structure (FFTW)
// planb      pointer to backward transform data structure (FFTW)
// ffttable   intermediate array used by the FFT routine (FFTPACK)
// ffttyp     type of FFT package; currently FFTPACK or FFTW

constexpr int maxprime = 15;
MDQC_EXTERN int iprime[3][maxprime];
MDQC_EXTERN int64_t planf;
MDQC_EXTERN int64_t planb;
MDQC_EXTERN std::vector<std::vector<double>> ffttable;
MDQC_EXTERN std::string ffttyp;
}
