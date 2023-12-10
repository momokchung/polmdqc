// Author: Moses KJ Chung
// Year:   2023

#include "mathUtils.h"
#include <cmath>
#include <vector>

namespace mathUtils{
//////////////////////////////////////////////////
//                                              //
//  symmetrize  --  symmetrize lower 2D matrix  //
//                                              //
//////////////////////////////////////////////////

void symmetrize(std::vector<std::vector<real>>& matrix)
{
    int N = matrix.size();
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            matrix[j][i] = matrix[i][j];
        }
    }
}

/////////////////////////////////////////////
//                                         //
//  dsyevd  --  wrapper for LAPACK dsyevd  //
//                                         //
/////////////////////////////////////////////

// Given matrix dimension N, row major 2d matrix A, and 1d matrix w,
// returns eigenvectors and eigenvalues of real symmetric matrix
extern "C" 
{
    extern void dsyevd_(char*, char*, int*, double*, int*, double*, double*, int*, int*, int*, int*);
}
void dsyevd(int n, real* A, real* w)
{
    char jobz = 'V';
    char uplo = 'U';
    int lda = n;
    int lwork = 1 + 6 * n + 2 * n * n;
    real* work = new real[lwork];
    int liwork = 3 + 5 * n;
    int* iwork = new int[liwork];
    int info;

    dsyevd_(&jobz, &uplo, &n, A, &lda, w, work, &lwork, iwork, &liwork, &info);

    delete[] work;
    delete[] iwork;
}

/////////////////////////////////////////
//                                     //
//  dgemm  --  wrapper for BLAS dgemm  //
//                                     //
/////////////////////////////////////////

// Given matrix A (M x K) and matrix B (K by N) computes C (M by N),
// where C = alpha * A * B + beta * C
extern "C"
{
    extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
}
void dgemm(int m, int n, int k, real* A, real* B, real* C, char transa, char transb)
{
    int lda = m;
    int ldb = k;
    if (transa == 'T')
        lda = k;
    if (transb == 'T')
        ldb = n;
    int ldc = m;
    real alpha = 1.;
    real beta = 0.;

    dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

}

///////////////////////////////////////////////////////////////
//                                                           //
//  triDiagSym  --  tridiagonalize real symmetric 2D matrix  //
//                                                           //
///////////////////////////////////////////////////////////////

// Burden R, Faires D, Burden A. Numerical Analysis 10th Ed. pg. 607-609
// Uses Householder algorithm to tridiagonalize real symmetric 2D matrix
void triDiagSym(std::vector<std::vector<real>>& A)
{
    // optimize A indexing
    int N = A.size();
    int Nm = N - 1;
    std::vector<real> v(N, 0.);
    std::vector<real> u(N, 0.);
    std::vector<real> z(N, 0.);

    // Step 1
    for (int k = 0; k < N - 2; ++k)
    {
        // Step 2
        real q = 0;
        for (int j = k + 1; j < N; ++j){
            q += A[j][k] * A[j][k];
        }

        // Step 3
        real alpha;
        if (fabs(A[k + 1][k]) < 1E-15)
        {
            alpha = -sqrt(q);
        }
        else
        {
            alpha = -sqrt(q) * A[k + 1][k] / fabs(A[k + 1][k]);
        }

        // Step 4
        real rsq = alpha * (alpha - A[k+1][k]);

        // Step 5
        v[k] = 0.;
        v[k + 1] = A[k + 1][k] - alpha;
        for (int j = k + 2; j < N; ++j)
        {
            v[j] = A[j][k];
        }

        // Step 6
        for (int j = k; j < N; ++j)
        {
            real utmp = 0;
            for (int i = k + 1; i < N; ++i)
            {
                utmp += A[j][i] * v[i];
            }
            u[j] = utmp / rsq;
        }

        // Step 7
        real prod = 0.;
        for (int i = k + 1; i < N; ++i)
        {
            prod += v[i] * u[i];
        }

        // Step 8
        for (int j = k; j < N; ++j)
        {
            z[j] = u[j] - prod / (2. * rsq) * v[j];
        }

        // Step 9
        for (int l = k + 1; l < N - 1; ++l)
        {
            // Step 10
            for (int j = l + 1; j < N; ++j)
            {
                A[j][l] = A[j][l] - v[l] * z[j] - v[j] * z[l];
                A[l][j] = A[j][l];
            }

            // Step 11
            A[l][l] = A[l][l] - 2. * v[l] * z[l];
        }
        
        // Step 12
        A[Nm][Nm] = A[Nm][Nm] - 2. * v[Nm] * z[Nm];

        // Step 13
        for (int j = k + 2; j < N; ++j)
        {
            A[k][j] = 0;
            A[j][k] = 0;
        }

        // Step 14
        A[k + 1][k] = A[k + 1][k] - v[k + 1] * z[k];
        A[k][k + 1] = A[k + 1][k];
    }
}
}
