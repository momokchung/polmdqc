//////////////////////////////////////////////////////
//                                                  //
//  tbFock.cpp  --  compute two body Fock elements  //
//                                                  //
//////////////////////////////////////////////////////


// #include "hartree.h"
// #include "katomsqm.h"
// #include "kbasis.h"
#include "kgbs.h"
// #include "kinetic.h"
// #include "kworker.h"
// #include "nuclear.h"
// #include "nuclearRepulsion.h"
// #include "overlap.h"
// #include "print.h"
// #include "valeev1.h"
#include <string>
#include <libint2.hpp>

namespace tbFock
{
void tbFock(std::vector<real>& G)
{
    std::string xyzfilename = "/Users/moseschung/polmdqc/example/QMexample/tmp.xyz";
    std::ifstream input_file(xyzfilename);
    std::vector<libint2::Atom> atoms = libint2::read_dotxyz(input_file);
    libint2::BasisSet obs(gbs::basisName, atoms);

    libint2::Engine engine(libint2::Operator::coulomb, obs.max_nprim(), obs.max_l());

    // libint2::Engine s_engine(libint2::Operator::overlap,  // will compute overlap ints
    //                         obs.max_nprim(),              // max # of primitives in shells this engine will accept
    //                         obs.max_l()                   // max angular momentum of shells this engine will accept
    //                         );
    // //
    auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
                                    // shell2bf[0] = index of the first basis function in shell 0
                                    // shell2bf[1] = index of the first basis function in shell 1
                                    // ...

    // buf[0] points to the target shell set after every call to engine.compute()
    const auto& buf = engine.results();
    // const auto& buf_vec = s_engine.results(); // will point to computed shell sets
    //                                           // const auto& is very important!

    // for(auto s1=0; s1!=obs.size(); ++s1) {
    //     for(auto s2=0; s2!=obs.size(); ++s2) {

    //         std::cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
    //         s_engine.compute(obs[s1], obs[s2]);
    //         std::cout << "done" << std::endl;
    //         auto ints_shellset = buf_vec[0];  // location of the computed integrals
    //         if (ints_shellset == nullptr)
    //             continue;  // nullptr returned if the entire shell-set was screened out

    //         auto bf1 = shell2bf[s1];  // first basis function in first shell
    //         auto n1 = obs[s1].size(); // number of basis functions in first shell
    //         auto bf2 = shell2bf[s2];  // first basis function in second shell
    //         auto n2 = obs[s2].size(); // number of basis functions in second shell

    //         // integrals are packed into ints_shellset in row-major (C) form
    //         // this iterates over integrals in this order
    //         for(auto f1=0; f1!=n1; ++f1)
    //             for(auto f2=0; f2!=n2; ++f2)
    //                 std::cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << std::endl;
    //     }
    // }
}
}
