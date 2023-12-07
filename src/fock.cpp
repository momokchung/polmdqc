///////////////////////////////////////////
//                                       //
//  fock.cpp  --  compute Fock elements  //
//                                       //
///////////////////////////////////////////


#include "katomsqm.h"
#include "fock.h"
#include "hartree.h"
#include "kbasis.h"
#include "kgbs.h"
#include "print.h"
#include <Eigen/Dense>
#include <string>

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

namespace fock
{
void fock_1body_libint(libint2::Operator obtype, std::vector<std::vector<real>>& mat)
{
    int N = basis::N;
    Matrix result(N, N);

    std::string xyzfilename = "/Users/moseschung/polmdqc/example/QCexample/tmp.xyz";
    std::ifstream input_file(xyzfilename);
    std::vector<libint2::Atom> atoms;
    for (int i = 0; i < atoms::n; ++i)
    {
        int an = atoms::core[i];
        real xi = atoms::coordx[i];
        real yi = atoms::coordy[i];
        real zi = atoms::coordz[i];
        libint2::Atom atomi{an,xi,yi,zi};
        atoms.push_back(atomi);
    }
    libint2::BasisSet obs(gbs::basisName, atoms);

    libint2::Engine engine(obtype, obs.max_nprim(), obs.max_l());

    if (obtype == libint2::Operator::nuclear)
        engine.set_params(make_point_charges(atoms));

    auto shell2bf = obs.shell2bf();
    const auto& buf = engine.results();

    for(auto s1=0; s1!=obs.size(); ++s1) {

        auto bf1 = shell2bf[s1]; // first basis function in this shell
        auto n1 = obs[s1].size();

        for(auto s2=0; s2<=s1; ++s2) {

            auto bf2 = shell2bf[s2];
            auto n2 = obs[s2].size();

            // compute shell pair
            engine.compute(obs[s1], obs[s2]);

            // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
            result.block(bf1, bf2, n1, n2) = buf_mat;
            if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
            result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mat[i][j] = result(i,j);
        }
    }
}

void fock_1body_lib()
{
    fock_1body_libint(libint2::Operator::overlap, hartree::S);
    fock_1body_libint(libint2::Operator::kinetic, hartree::KE);
    fock_1body_libint(libint2::Operator::nuclear, hartree::NE);
}

void fock_2body_lib()
{
    int N = basis::N;

    Matrix G = Matrix::Zero(N,N);
    Matrix D = Matrix::Zero(N,N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            D(i,j) = hartree::D[N*i+j];
        }
    }

    std::string xyzfilename = "/Users/moseschung/polmdqc/example/QCexample/tmp.xyz";
    std::ifstream input_file(xyzfilename);
    std::vector<libint2::Atom> atoms;
    for (int i = 0; i < atoms::n; ++i)
    {
        int an = atoms::core[i];
        real xi = atoms::coordx[i];
        real yi = atoms::coordy[i];
        real zi = atoms::coordz[i];
        libint2::Atom atomi{an,xi,yi,zi};
        atoms.push_back(atomi);
    }
    libint2::BasisSet obs(gbs::basisName, atoms);

    libint2::Engine engine(libint2::Operator::coulomb, obs.max_nprim(), obs.max_l());

    auto shell2bf = obs.shell2bf();

    // buf[0] points to the target shell set after every call to engine.compute()
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for(auto s1=0; s1!=obs.size(); ++s1) {

        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1 = obs[s1].size();   // number of basis functions in this shell

        for(auto s2=0; s2<=s1; ++s2) {

            auto bf2_first = shell2bf[s2];
            auto n2 = obs[s2].size();

            for(auto s3=0; s3<=s1; ++s3) {

                auto bf3_first = shell2bf[s3];
                auto n3 = obs[s3].size();

                const auto s4_max = (s1 == s3) ? s2 : s3;
                for(auto s4=0; s4<=s4_max; ++s4) {

                    auto bf4_first = shell2bf[s4];
                    auto n4 = obs[s4].size();

                    // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                    auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                    auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                    auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                    auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                    engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
                    const auto* buf_1234 = buf[0];
                    if (buf_1234 == nullptr)
                        continue; // if all integrals screened out, skip to next quartet

                    // ANSWER
                    // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
                    //    F(a,b) += (ab|cd) * D(c,d)
                    //    F(c,d) += (ab|cd) * D(a,b)
                    //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
                    //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
                    //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
                    //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
                    // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
                    //    i.e. the number of the integrals/sets equivalent to it
                    // 3) the end result must be symmetrized
                    for(auto f1=0, f1234=0; f1!=n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for(auto f2=0; f2!=n2; ++f2) {
                            const auto bf2 = f2 + bf2_first;
                            for(auto f3=0; f3!=n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                                    const auto bf4 = f4 + bf4_first;

                                    const auto value = buf_1234[f1234];

                                    const auto value_scal_by_deg = value * s1234_deg;

                                    G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                                    G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                                    G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                                    G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                                    G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                                    G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // symmetrize the result and return
    Matrix Gt = G.transpose();
    G = 0.5 * (G + Gt);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            hartree::G[N*i+j] = G(i,j);
        }
    }
}
}
