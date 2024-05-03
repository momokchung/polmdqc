#pragma once
#include "precision.h"
#include <vector>

namespace polmdqc
{
namespace spacefill1
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;

    extern std::vector<real> surf;
    extern std::vector<real> vol;
    extern std::vector<real> mean;
    extern std::vector<real> gauss;
    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill2
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;

    extern std::vector<real> surf;
    extern std::vector<real> vol;
    extern std::vector<real> mean;
    extern std::vector<real> gauss;
    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill3
{
    extern real eps;
    extern real eps2;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;

    extern std::vector<real> surf;
    extern std::vector<real> vol;
    extern std::vector<real> mean;
    extern std::vector<real> gauss;
    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill4
{
    extern real eps;
    extern real eps2;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;

    extern std::vector<real> surf;
    extern std::vector<real> vol;
    extern std::vector<real> mean;
    extern std::vector<real> gauss;
    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill5
{
    extern real eps;
    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill6
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;
    extern real wsurf2;
    extern real wvol2;
    extern real wmean2;
    extern real wgauss2;
    extern real wsurf3;
    extern real wvol3;
    extern real wmean3;
    extern real wgauss3;
    extern real wsurf4;
    extern real wvol4;
    extern real wmean4;
    extern real wgauss4;
}

namespace spacefill7
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;
    extern real wsurf2;
    extern real wvol2;
    extern real wmean2;
    extern real wgauss2;
    extern real wsurf3;
    extern real wvol3;
    extern real wmean3;
    extern real wgauss3;
    extern real wsurf4;
    extern real wvol4;
    extern real wmean4;
    extern real wgauss4;
}

namespace spacefill8
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;
    extern real wsurf2;
    extern real wvol2;
    extern real wmean2;
    extern real wgauss2;
    extern real wsurf3;
    extern real wvol3;
    extern real wmean3;
    extern real wgauss3;
    extern real wsurf4;
    extern real wvol4;
    extern real wmean4;
    extern real wgauss4;
}

namespace spacefill9
{
    extern real eps;
    extern std::vector<real> surf;
    extern std::vector<real> vol;
    extern std::vector<real> mean;
    extern std::vector<real> gauss;
}

namespace spacefill10
{
    extern real eps;
    extern real epsg;
    extern real eps2;

    extern std::vector<real> dsurf;
    extern std::vector<real> dvol;
    extern std::vector<real> dmean;
    extern std::vector<real> dgauss;
}

namespace spacefill11
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;
    extern real wsurf2;
    extern real wvol2;
    extern real wmean2;
    extern real wgauss2;
}

namespace spacefill12
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;
}

namespace spacefill13
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;
}

namespace spacefill14
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
}

namespace spacefill15
{
    extern real eps;
    extern real wsurf;
    extern real wvol;
    extern real wmean;
    extern real wgauss;
}

namespace spacefill16
{
    extern real eps;
    extern real wsurf0;
    extern real wvol0;
    extern real wmean0;
    extern real wgauss0;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;

    // // current AlphaMol can't handle linear structure
    // extern real wsurf2;
    // extern real wvol2;
    // extern real wmean2;
    // extern real wgauss2;

    // // current AlphaMol can't handle linear structure
    // extern real wsurf3;
    // extern real wvol3;
    // extern real wmean3;
    // extern real wgauss3;

    // // current AlphaMol can't handle linear structure
    // extern real wsurf4;
    // extern real wvol4;
    // extern real wmean4;
    // extern real wgauss4;

    // // current AlphaMol can't handle linear structure
    // extern real wsurf5;
    // extern real wvol5;
    // extern real wmean5;
    // extern real wgauss5;
}

namespace spacefill17
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;

    // // current AlphaMol can't handle planar cases
    // extern real wsurf2;
    // extern real wvol2;
    // extern real wmean2;
    // extern real wgauss2;

    // // current AlphaMol can't handle planar cases
    // extern real wsurf3;
    // extern real wvol3;
    // extern real wmean3;
    // extern real wgauss3;
}

namespace spacefill18
{
    extern real eps;
    extern real wsurf1;
    extern real wvol1;
    extern real wmean1;
    extern real wgauss1;
    extern real wsurf2;
    extern real wvol2;
    extern real wmean2;
    extern real wgauss2;
    extern real wsurf3;
    extern real wvol3;
    extern real wmean3;
    extern real wgauss3;

    // // current AlphaMol can't near symmetric cases
    // extern real wsurf4;
    // extern real wvol4;
    // extern real wmean4;
    // extern real wgauss4;
}
}