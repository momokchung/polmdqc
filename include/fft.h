////////////////////////////////////////////////////////
//                                                    //
//  fft.h  --  Fast Fourier transform control values  //
//                                                    //
////////////////////////////////////////////////////////

// maxprime   maximum number of prime factors of FFT dimension
// iprime     prime factorization of each FFT dimension (FFTPACK)
// planf      pointer to forward transform data structure (FFTW)
// planb      pointer to backward transform data structure (FFTW)
// ffttable   intermediate array used by the FFT routine (FFTPACK)
// ffttyp     type of FFT package; currently FFTPACK or FFTW


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN constexpr int maxprime = 15;
QCMD_EXTERN int iprime[3][maxprime];
QCMD_EXTERN int64_t planf;
QCMD_EXTERN int64_t planb;
QCMD_EXTERN std::vector<std::vector<double>> ffttable;
QCMD_EXTERN std::string ffttyp;
