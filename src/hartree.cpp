/////////////////////////////////////////////
//                                         //
//  hartree.cpp  --  compute Hartree Fock  //
//                                         //
/////////////////////////////////////////////


#include "hartree.h"
#include "katomsqm.h"
#include "kbasis.h"
#include "kgbs.h"
#include "kinetic.h"
#include "kworker.h"
#include "nuclear.h"
#include "nuclearRepulsion.h"
#include "overlap.h"
#include <iostream>
#include <libint2.hpp>

namespace hartree
{
// s_tolerance    tolerance for symmetric vs canonical orthogonalization
// s_i            starting index because of tolerance

real s_tolerance = 1E-7;
int s_i;

// H     column major Hamiltonian matrix
// F     column major Fock matrix
// Fp    transformed F
// C     coefficient matrix
// E     energy eigenvalues of FC = CE
// D     density matrix
// X     column major S^(-1/2) matrix
// vS    column major eigenvector matrix of S
// wS    eigenvalues of S

std::vector<real> H;
std::vector<real> F;
std::vector<real> Fp;
std::vector<real> C;
std::vector<real> E;
std::vector<real> D;
std::vector<real> X;
std::vector<real> vS;
std::vector<real> wS;


/////////////////////////////////////////////////////
//                                                 //
//  void rhf  --  compute restricted Hartree Fock  //
//                                                 //
/////////////////////////////////////////////////////

// current implementation configuration: scf_incore, convergence_conventional, guess_core
// need to implement scf_direc, scf_conventional, scf_densityFit, convergence_DIIS, guess_sad
// need to implement other convergence methods as well (damping, 2nd order, etc.)
void rhf()
{
    // libint2::initialize();
    int N = basis::N;
    int N2 = N * N;
    int nElec = atoms::nElec;
    // check if number of electrons are even
    if ( nElec % 2 != 0)
    {
        std::cerr << "Number of electrons is not even, cannot do RHF " << "\n";
        std::exit(1);
    }
    int ndocc = nElec/2;

    // allocate and initialize matrices
    H.resize(0, 0.);
    F.resize(0);
    Fp.resize(0);
    C.resize(0);
    E.resize(0);
    D.resize(0);
    H.resize(N2, 0.);
    F.reserve(N2);
    Fp.reserve(N2);
    C.reserve(N2);
    E.reserve(N);
    D.reserve(N2);

    // construct overlap matrix
    // change to overlap() later
    overlap::overlapOS();
    std::vector<std::vector<real>>& S = (gbs::basisType == gbs::BasisType::cartesian) ? overlap::cartS : overlap::sphS;

    // diagonalize overlap matrix
    eigenS(S);

    // construct KE matrix
    // change to kinetic() later
    kinetic::kineticOS();
    std::vector<std::vector<real>>& KE = (gbs::basisType == gbs::BasisType::cartesian) ? kinetic::cartKE : kinetic::sphKE;

    // construct NE matrix
    // change to nuclear() later
    nuclear::nuclearOS();
    std::vector<std::vector<real>>& NE = (gbs::basisType == gbs::BasisType::cartesian) ? nuclear::cartNE : nuclear::sphNE;

    // compute nuclear nuclear repulsion
    nuclearRepulsion::nuclearRepulsion();

    // form the Hamiltonian matrix
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            H[N * i + j] = KE[i][j] + NE[i][j];
        }
    }

    // core guess
    guess();

    // conventional scf
    // libint2::finalize();
}


////////////////////////////////////////////////////////////////////
//                                                                //
//  void eigenS  --  solve the eigenvalues and eigenvectors of S  //
//                                                                //
////////////////////////////////////////////////////////////////////

// do symmetric orthogonalization if all eigenvalues are smaller than s_tolerance
// else do canonical orthogonalization
void eigenS(std::vector<std::vector<real>>& S)
{
    int N = basis::N;

    // allocate memory for vS and wS
    vS.resize(0);
    wS.resize(0);
    vS.reserve(N * N);
    wS.reserve(N);

    // flatten S to vS
    mathUtils::flatten(N, N, S, vS);

    // solve eigenvalues and eigenvectors of S
    mathUtils::dsyevd(N, vS.data(), wS.data());

    s_i = 0;
    while (wS[s_i] < s_tolerance)
    {
        s_i += 1;
    }

    if (s_i == 0) // symmetric orthogonalization
    {
        // allocate and initialize X
        X.resize(N * N, 0.);

        // assign diagonal elements of X
        for (int i = 0; i < N; ++i)
        {
            X[N * i + i] = 1/sqrt(wS[i]);
        }

        // reference workerN2
        std::vector<real>& workerN2 = worker::workerN2_1;
        
        // now we just need to do U * X * U^T!
        // U * X
        mathUtils::dgemm(N, N, N, vS.data(), X.data(), workerN2.data());
        // (U * X) * U^T
        mathUtils::dgemm(N, N, N, workerN2.data(), vS.data(), X.data(), 'N', 'T');
    }
    else // canonical orthogonalization
    {
        // allocate and initialize X
        int Nmsi = N - s_i;
        X.resize(Nmsi * N, 0.);
        for (int i = 0; i < Nmsi; ++i)
        {
            real s12 = sqrt(wS[i + s_i]);
            for (int j = 0; j < N; ++j)
            {
                X[N * i + j] = vS[N * (i + s_i) + j] / s12;
            }
        }
    }

    // print to debug
    // printf("[");
    // for (int i = 0; i < N; ++i)
    // {
    //     printf("%10.5f,", wS[i]);
    // }
    // printf("]\n");
    // for (int i = 0; i < N; ++i)
    // {
    //     printf("[");
    //     for (int j = 0; j < (N - s_i); ++j)
    //     {
    //         printf("%20.16f,", X[i + j * N]);
    //     }
    //     printf("],\n");
    // }
}


////////////////////////////////////////////////////////
//                                                    //
//  void guess  --  guess the initial density matrix  //
//                                                    //
////////////////////////////////////////////////////////

// current implementation is core, need to implement SAD in the future
void guess()
{
    int N = basis::N;
    int nElec = atoms::nElec;
    int ndocc = nElec/2;

    // dimension of X is N x (N - S_i)
    int Xn = N - s_i;

    // reference workerN2
    std::vector<real>& workerN2 = worker::workerN2_1;

    // transformed Fock matrix F' = X^T * H * X
    // X^T * H
    mathUtils::dgemm(Xn, N, N, X.data(), H.data(), workerN2.data(), 'T', 'N');
    // (X^T * H) * X
    mathUtils::dgemm(Xn, Xn, N, workerN2.data(), X.data(), Fp.data(), 'N', 'N');

    // solve eigenvalues and eigenvectors of Fp, Fp now becomes Cp
    mathUtils::dsyevd(Xn, Fp.data(), E.data());

    // transform Cp back to C, C = X * Cp
    mathUtils::dgemm(N, Xn, Xn, X.data(), Fp.data(), C.data());
    
    // build density matrix D = C * C^T
    mathUtils::dgemm(N, N, ndocc, C.data(), C.data(), D.data(), 'N', 'T');

    // mathUtils::dgemm(N, N, N, D.data(), H.data(), workerN2.data(), 'N', 'N');
    // real energy = 0;
    // for (int i = 0; i < N; ++i)
    // {
    //     energy += workerN2[N * i + i];
    // }
    // printf("%18.14f\n", energy);
}


///////////////////////////////////////////
//                                       //
//  void scf  --  perform scf procedure  //
//                                       //
///////////////////////////////////////////

// current implementation is scf_incore + convergence_conventional,
// need to implement scf_direct, scf_conventional, and scf_densityFit
// need to implement convergence_DIIS, convergence_damping, and other convergence methods
void scf()
{
    
}
}