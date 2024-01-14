#pragma once
#include "precision.h"
#include <vector>

namespace polmdqc
{
namespace hippo1
{
    real eps = 1e-10;

    int nem = 9;

    real einter = -7.3522485663420181E+00;
    real esum = -7.3522485663420181E+00;
    real em = -7.3522485663420181E+00;

    std::vector<real> aesum = {
        1.3692680424181254E+00,
       -7.9652177693810611E-01,
       -4.2488705486510279E+00,
       -4.5558705761474814E+00,
        4.3991005456822418E-01,
        4.3983623840824881E-01,
    };
    std::vector<real> aem = {
        1.3692680424181254E+00,
       -7.9652177693810611E-01,
       -4.2488705486510279E+00,
       -4.5558705761474814E+00,
        4.3991005456822418E-01,
        4.3983623840824881E-01,
    };
}

namespace hippo2
{
    real eps1 = 1e-10;
    real eps2 = 1e-6;

    std::vector<std::vector<real>> desum = {
        {-8.6567440108873439E-01,  7.9536557077020298E-01, -1.6681990701249587E-02,},
        {-1.2239035280358608E+00,  1.8430306353002854E-01, -1.3547115567577107E-03,},
        {-1.2471016517145934E+01, -3.5660218891848977E-01,  3.7766712757924870E-02,},
        { 1.4153204919015218E+01, -3.1242203751583215E+00,  4.4163231284820303E-02,},
        { 2.0420385795462184E-01,  1.2430178733007831E+00, -3.4849988956922684E-01,},
        { 2.0318566930068915E-01,  1.2581360564757966E+00,  2.8460664778448891E-01,},
    };

    std::vector<std::vector<real>> dem = {
        {-8.6567440108873439E-01,  7.9536557077020298E-01, -1.6681990701249587E-02,},
        {-1.2239035280358608E+00,  1.8430306353002854E-01, -1.3547115567577107E-03,},
        {-1.2471016517145934E+01, -3.5660218891848977E-01,  3.7766712757924870E-02,},
        { 1.4153204919015218E+01, -3.1242203751583215E+00,  4.4163231284820303E-02,},
        { 2.0420385795462184E-01,  1.2430178733007831E+00, -3.4849988956922684E-01,},
        { 2.0318566930068915E-01,  1.2581360564757966E+00,  2.8460664778448891E-01,},
    };
}
}
