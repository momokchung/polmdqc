#include "alfp.h"
#include "alphmol.h"
#include "final.h"
#include "inform.h"
#include "spacefill.h"
#include "spacefillT.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("spacefill-1", "[spacefill][AMOEBA][water09]") {
    // Van der Waals Area and Volume
    // include hydrogen
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/water09.key",
        "../../test/testFiles/spacefill/water09.xyz",
        "1","Y","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    REQUIRE((alfdigit % 2) == 0);

    COMPARE_REALS(wsurf, spacefill1::wsurf, spacefill1::eps);
    COMPARE_REALS(wvol, spacefill1::wvol, spacefill1::eps);
    COMPARE_REALS(wmean, spacefill1::wmean, spacefill1::eps);
    COMPARE_REALS(wgauss, spacefill1::wgauss, spacefill1::eps);
    
    int n = 6;
    for (int i = 0; i < n; i++) {
        COMPARE_REALS(surf[i], spacefill1::surf[i], spacefill1::eps);
        COMPARE_REALS(vol[i], spacefill1::vol[i], spacefill1::eps);
        COMPARE_REALS(mean[i], spacefill1::mean[i], spacefill1::eps);
        COMPARE_REALS(gauss[i], spacefill1::gauss[i], spacefill1::eps);
    }

    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill1::dsurf[i], spacefill1::eps);
        COMPARE_REALS(dvol[i], spacefill1::dvol[i], spacefill1::eps);
        COMPARE_REALS(dmean[i], spacefill1::dmean[i], spacefill1::eps);
        COMPARE_REALS(dgauss[i], spacefill1::dgauss[i], spacefill1::eps);
    }

    final();
}

TEST_CASE("spacefill-2", "[spacefill][AMOEBA][water09]") {
    // Van der Waals Area and Volume
    // exclude hydrogen
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/water09.key",
        "../../test/testFiles/spacefill/water09.xyz",
        "1","N","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill2::wsurf, spacefill2::eps);
    COMPARE_REALS(wvol, spacefill2::wvol, spacefill2::eps);
    COMPARE_REALS(wmean, spacefill2::wmean, spacefill2::eps);
    COMPARE_REALS(wgauss, spacefill2::wgauss, spacefill2::eps);

    int n = 6;
    for (int i = 0; i < n; i++) {
        COMPARE_REALS(surf[i], spacefill2::surf[i], spacefill2::eps);
        COMPARE_REALS(vol[i], spacefill2::vol[i], spacefill2::eps);
        COMPARE_REALS(mean[i], spacefill2::mean[i], spacefill2::eps);
        COMPARE_REALS(gauss[i], spacefill2::gauss[i], spacefill2::eps);
    }

    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill2::dsurf[i], spacefill2::eps);
        COMPARE_REALS(dvol[i], spacefill2::dvol[i], spacefill2::eps);
        COMPARE_REALS(dmean[i], spacefill2::dmean[i], spacefill2::eps);
        COMPARE_REALS(dgauss[i], spacefill2::dgauss[i], spacefill2::eps);
    }

    final();
}

TEST_CASE("spacefill-3", "[spacefill][AMOEBA][water09]") {
    // Accessible Area and Excluded Volume
    // include hydrogen
    int argc = 10;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/water09.key",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","Y","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill3::wsurf, spacefill3::eps);
    COMPARE_REALS(wvol, spacefill3::wvol, spacefill3::eps);
    COMPARE_REALS(wmean, spacefill3::wmean, spacefill3::eps);
    COMPARE_REALS(wgauss, spacefill3::wgauss, spacefill3::eps);

    int n = 6;
    for (int i = 0; i < n; i++) {
        COMPARE_REALS(surf[i], spacefill3::surf[i], spacefill3::eps);
        COMPARE_REALS(vol[i], spacefill3::vol[i], spacefill3::eps);
        COMPARE_REALS(mean[i], spacefill3::mean[i], spacefill3::eps);
        COMPARE_REALS(gauss[i], spacefill3::gauss[i], spacefill3::eps);
    }

    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill3::dsurf[i], spacefill3::eps2);
        COMPARE_REALS(dvol[i], spacefill3::dvol[i], spacefill3::eps2);
        COMPARE_REALS(dmean[i], spacefill3::dmean[i], spacefill3::eps2);
        COMPARE_REALS(dgauss[i], spacefill3::dgauss[i], spacefill3::eps2);
    }

    final();
}

TEST_CASE("spacefill-4", "[spacefill][AMOEBA][water09]") {
    // Accessible Area and Excluded Volume
    // exclude hydrogen
    int argc = 10;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/water09.key",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","N","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill4::wsurf, spacefill4::eps);
    COMPARE_REALS(wvol, spacefill4::wvol, spacefill4::eps);
    COMPARE_REALS(wmean, spacefill4::wmean, spacefill4::eps);
    COMPARE_REALS(wgauss, spacefill4::wgauss, spacefill4::eps);

    int n = 6;
    for (int i = 0; i < n; i++) {
        COMPARE_REALS(surf[i], spacefill4::surf[i], spacefill4::eps);
        COMPARE_REALS(vol[i], spacefill4::vol[i], spacefill4::eps);
        COMPARE_REALS(mean[i], spacefill4::mean[i], spacefill4::eps);
        COMPARE_REALS(gauss[i], spacefill4::gauss[i], spacefill4::eps);
    }

    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill4::dsurf[i], spacefill4::eps2);
        COMPARE_REALS(dvol[i], spacefill4::dvol[i], spacefill4::eps2);
        COMPARE_REALS(dmean[i], spacefill4::dmean[i], spacefill4::eps2);
        COMPARE_REALS(dgauss[i], spacefill4::dgauss[i], spacefill4::eps2);
    }

    final();
}

TEST_CASE("spacefill-5", "[spacefill][AMOEBA][water09]") {
    // numerical gradient
    int argc = 11;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/water09.key",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","Y","N","N","Y","1e-5"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    int n = 6;
    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(ndsurf[i], spacefill5::dsurf[i], spacefill5::eps);
        COMPARE_REALS(ndvol[i], spacefill5::dvol[i], spacefill5::eps);
        COMPARE_REALS(ndmean[i], spacefill5::dmean[i], spacefill5::eps);
        COMPARE_REALS(ndgauss[i], spacefill5::dgauss[i], spacefill5::eps);
    }

    final();
}

TEST_CASE("spacefill-6", "[spacefill][AMOEBA][alatet]") {
    // 1
    int argc1 = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 9;
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 10;
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 10;
    const char* strings4[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill6::wsurf1, spacefill6::eps);
    COMPARE_REALS(wvol, spacefill6::wvol1, spacefill6::eps);
    COMPARE_REALS(wmean, spacefill6::wmean1, spacefill6::eps);
    COMPARE_REALS(wgauss, spacefill6::wgauss1, spacefill6::eps);
    final();

    spacefill(argc2, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill6::wsurf2, spacefill6::eps);
    COMPARE_REALS(wvol, spacefill6::wvol2, spacefill6::eps);
    COMPARE_REALS(wmean, spacefill6::wmean2, spacefill6::eps);
    COMPARE_REALS(wgauss, spacefill6::wgauss2, spacefill6::eps);
    final();

    spacefill(argc3, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill6::wsurf3, spacefill6::eps);
    COMPARE_REALS(wvol, spacefill6::wvol3, spacefill6::eps);
    COMPARE_REALS(wmean, spacefill6::wmean3, spacefill6::eps);
    COMPARE_REALS(wgauss, spacefill6::wgauss3, spacefill6::eps);
    final();

    spacefill(argc4, argv4);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill6::wsurf4, spacefill6::eps);
    COMPARE_REALS(wvol, spacefill6::wvol4, spacefill6::eps);
    COMPARE_REALS(wmean, spacefill6::wmean4, spacefill6::eps);
    COMPARE_REALS(wgauss, spacefill6::wgauss4, spacefill6::eps);
    final();
}

TEST_CASE("spacefill-7", "[spacefill][AMOEBA][alatet]") {
    // inactive
    // 1
    int argc1 = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 9;
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 10;
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 10;
    const char* strings4[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill7::wsurf1, spacefill7::eps);
    COMPARE_REALS(wvol, spacefill7::wvol1, spacefill7::eps);
    COMPARE_REALS(wmean, spacefill7::wmean1, spacefill7::eps);
    COMPARE_REALS(wgauss, spacefill7::wgauss1, spacefill7::eps);
    final();

    spacefill(argc2, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill7::wsurf2, spacefill7::eps);
    COMPARE_REALS(wvol, spacefill7::wvol2, spacefill7::eps);
    COMPARE_REALS(wmean, spacefill7::wmean2, spacefill7::eps);
    COMPARE_REALS(wgauss, spacefill7::wgauss2, spacefill7::eps);
    final();

    spacefill(argc3, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill7::wsurf3, spacefill7::eps);
    COMPARE_REALS(wvol, spacefill7::wvol3, spacefill7::eps);
    COMPARE_REALS(wmean, spacefill7::wmean3, spacefill7::eps);
    COMPARE_REALS(wgauss, spacefill7::wgauss3, spacefill7::eps);
    final();

    spacefill(argc4, argv4);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill7::wsurf4, spacefill7::eps);
    COMPARE_REALS(wvol, spacefill7::wvol4, spacefill7::eps);
    COMPARE_REALS(wmean, spacefill7::wmean4, spacefill7::eps);
    COMPARE_REALS(wgauss, spacefill7::wgauss4, spacefill7::eps);
    final();
}

TEST_CASE("spacefill-8", "[spacefill][AMOEBA][alatet]") {
    // active
    // 1
    int argc1 = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 9;
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 10;
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 10;
    const char* strings4[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill8::wsurf1, spacefill8::eps);
    COMPARE_REALS(wvol, spacefill8::wvol1, spacefill8::eps);
    COMPARE_REALS(wmean, spacefill8::wmean1, spacefill8::eps);
    COMPARE_REALS(wgauss, spacefill8::wgauss1, spacefill8::eps);
    final();

    spacefill(argc2, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill8::wsurf2, spacefill8::eps);
    COMPARE_REALS(wvol, spacefill8::wvol2, spacefill8::eps);
    COMPARE_REALS(wmean, spacefill8::wmean2, spacefill8::eps);
    COMPARE_REALS(wgauss, spacefill8::wgauss2, spacefill8::eps);
    final();

    spacefill(argc3, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill8::wsurf3, spacefill8::eps);
    COMPARE_REALS(wvol, spacefill8::wvol3, spacefill8::eps);
    COMPARE_REALS(wmean, spacefill8::wmean3, spacefill8::eps);
    COMPARE_REALS(wgauss, spacefill8::wgauss3, spacefill8::eps);
    final();

    spacefill(argc4, argv4);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill8::wsurf4, spacefill8::eps);
    COMPARE_REALS(wvol, spacefill8::wvol4, spacefill8::eps);
    COMPARE_REALS(wmean, spacefill8::wmean4, spacefill8::eps);
    COMPARE_REALS(wgauss, spacefill8::wgauss4, spacefill8::eps);
    final();
}

TEST_CASE("spacefill-9", "[spacefill][AMOEBA][alatet]") {
    // atomic surface area, volume, mean curvature, and gaussian curvature
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    
    int n = 42;
    for (int i = 0; i < n; i++) {
        COMPARE_REALS(surf[i], spacefill9::surf[i], spacefill9::eps);
        COMPARE_REALS(vol[i], spacefill9::vol[i], spacefill9::eps);
        COMPARE_REALS(mean[i], spacefill9::mean[i], spacefill9::eps);
        COMPARE_REALS(gauss[i], spacefill9::gauss[i], spacefill9::eps);
    }

    final();
}

TEST_CASE("spacefill-10", "[spacefill][AMOEBA][alatet]") {
    // gradient surface area, volume, mean curvature, and gaussian curvature
    int argc = 10;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatet.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","Y","Y","1e-5"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    
    int n = 42;
    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill10::dsurf[i], spacefill10::eps);
        COMPARE_REALS(dvol[i], spacefill10::dvol[i], spacefill10::eps);
        COMPARE_REALS(dmean[i], spacefill10::dmean[i], spacefill10::eps);
        COMPARE_REALS(dgauss[i], spacefill10::dgauss[i], spacefill10::epsg);
    }

    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(ndsurf[i], spacefill10::dsurf[i], spacefill10::eps2);
        COMPARE_REALS(ndvol[i], spacefill10::dvol[i], spacefill10::eps2);
        COMPARE_REALS(ndmean[i], spacefill10::dmean[i], spacefill10::eps2);
        COMPARE_REALS(ndgauss[i], spacefill10::dgauss[i], spacefill10::eps2);
    }

    final();
}

TEST_CASE("spacefill-11", "[spacefill][AMOEBA][concat]") {
    // ability to read concatenated files
    int argc = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/alatetwater.key",
        "../../test/testFiles/spacefill/alatetwater.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/wateralatet.key",
        "../../test/testFiles/spacefill/wateralatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);

    test = true;

    spacefill(argc, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill11::wsurf1, spacefill11::eps);
    COMPARE_REALS(wvol, spacefill11::wvol1, spacefill11::eps);
    COMPARE_REALS(wmean, spacefill11::wmean1, spacefill11::eps);
    COMPARE_REALS(wgauss, spacefill11::wgauss1, spacefill11::eps);
    final();

    spacefill(argc, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill11::wsurf2, spacefill11::eps);
    COMPARE_REALS(wvol, spacefill11::wvol2, spacefill11::eps);
    COMPARE_REALS(wmean, spacefill11::wmean2, spacefill11::eps);
    COMPARE_REALS(wgauss, spacefill11::wgauss2, spacefill11::eps);
    final();
}

TEST_CASE("spacefill-12", "[spacefill][AMOEBA][chloride27]") {
    // symmetric object
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/chloride27.key",
        "../../test/testFiles/spacefill/chloride27.xyz",
        "1","Y","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);
    int argc2 = 9;
    const char *strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/chloride27_2.key",
        "../../test/testFiles/spacefill/chloride27_2.xyz",
        "1","Y" "N","N","N"};
    char** argv2 = const_cast<char**>(strings2);

    test = true;

    spacefill(argc, argv);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill12::wsurf, spacefill12::eps);
    COMPARE_REALS(wvol, spacefill12::wvol, spacefill12::eps);
    COMPARE_REALS(wmean, spacefill12::wmean, spacefill12::eps);
    COMPARE_REALS(wgauss, spacefill12::wgauss, spacefill12::eps);
    int n = 27;
    for (int i = 0; i < 3*n; i++) {
        COMPARE_REALS(dsurf[i], spacefill12::dsurf[i], spacefill12::eps2);
        COMPARE_REALS(dvol[i], spacefill12::dvol[i], spacefill12::eps2);
    }
    final();

    spacefill(argc2, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill12::wsurf2, spacefill12::eps2);
    COMPARE_REALS(wvol, spacefill12::wvol2, spacefill12::eps2);
    COMPARE_REALS(wmean, spacefill12::wmean2, spacefill12::eps2);
    COMPARE_REALS(wgauss, spacefill12::wgauss2, spacefill12::eps2);
    final();
}

TEST_CASE("spacefill-13", "[spacefill][AMOEBA][3ibk]") {
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/3ibk.key",
        "../../test/testFiles/spacefill/3ibk.xyz",
        "1","Y","N","N","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill13::wsurf, spacefill13::eps);
    COMPARE_REALS(wvol, spacefill13::wvol, spacefill13::eps);
    COMPARE_REALS(wmean, spacefill13::wmean, spacefill13::eps);
    COMPARE_REALS(wgauss, spacefill13::wgauss, spacefill13::eps);

    final();
}

TEST_CASE("spacefill-14", "[spacefill][AMOEBA][3cln]") {
    int argc = 10;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/3cln.key",
        "../../test/testFiles/spacefill/3cln.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill14::wsurf, spacefill14::eps);
    COMPARE_REALS(wvol, spacefill14::wvol, spacefill14::eps);

    final();
}

TEST_CASE("spacefill-15", "[spacefill][AMOEBA][waterbox30]") {
    int argc = 9;
    const char* strings[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/waterbox30.key",
        "../../test/testFiles/spacefill/waterbox30.xyz",
        "1","Y","N","N","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    spacefill(argc, argv);

    REQUIRE(alfmeth == AlfMethod::AlphaMol);

    COMPARE_REALS(wsurf, spacefill15::wsurf, spacefill15::eps);
    COMPARE_REALS(wvol, spacefill15::wvol, spacefill15::eps);
    COMPARE_REALS(wmean, spacefill15::wmean, spacefill15::eps);
    COMPARE_REALS(wgauss, spacefill15::wgauss, spacefill15::eps);

    final();
}

TEST_CASE("spacefill-16", "[spacefill][AMOEBA][lchloride]") {
    int argc = 9;
    const char* strings0[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/chloride1.key",
        "../../test/testFiles/spacefill/chloride1.xyz",
        "1","Y","N","N","N"
    };
    char** argv0 = const_cast<char**>(strings0);
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/lchloride2.key",
        "../../test/testFiles/spacefill/lchloride2.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/lchloride3.key",
        "../../test/testFiles/spacefill/lchloride3.xyz",
        "1","Y","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/lchloride4.key",
        "../../test/testFiles/spacefill/lchloride4.xyz",
        "1","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    const char* strings4[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/lchloride5.key",
        "../../test/testFiles/spacefill/lchloride5.xyz",
        "1","Y","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);
    const char* strings5[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/lchloride6.key",
        "../../test/testFiles/spacefill/lchloride6.xyz",
        "1","Y","N","N","N"
    };
    char** argv5 = const_cast<char**>(strings5);

    test = true;

    spacefill(argc, argv0);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf0, spacefill16::eps);
    COMPARE_REALS(wvol, spacefill16::wvol0, spacefill16::eps);
    COMPARE_REALS(wmean, spacefill16::wmean0, spacefill16::eps);
    COMPARE_REALS(wgauss, spacefill16::wgauss0, spacefill16::eps);
    final();

    spacefill(argc, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf1, spacefill16::eps);
    COMPARE_REALS(wvol, spacefill16::wvol1, spacefill16::eps);
    COMPARE_REALS(wmean, spacefill16::wmean1, spacefill16::eps);
    COMPARE_REALS(wgauss, spacefill16::wgauss1, spacefill16::eps);
    final();

    spacefill(argc, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf2, spacefill16::eps2);
    COMPARE_REALS(wvol, spacefill16::wvol2, spacefill16::eps2);
    COMPARE_REALS(wmean, spacefill16::wmean2, spacefill16::eps2);
    COMPARE_REALS(wgauss, spacefill16::wgauss2, spacefill16::eps2);
    final();

    spacefill(argc, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf3, spacefill16::eps2);
    COMPARE_REALS(wvol, spacefill16::wvol3, spacefill16::eps2);
    COMPARE_REALS(wmean, spacefill16::wmean3, spacefill16::eps2);
    COMPARE_REALS(wgauss, spacefill16::wgauss3, spacefill16::eps2);
    final();

    spacefill(argc, argv4);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf4, spacefill16::eps2);
    COMPARE_REALS(wvol, spacefill16::wvol4, spacefill16::eps2);
    COMPARE_REALS(wmean, spacefill16::wmean4, spacefill16::eps2);
    COMPARE_REALS(wgauss, spacefill16::wgauss4, spacefill16::eps2);
    final();

    spacefill(argc, argv5);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill16::wsurf5, spacefill16::eps2);
    COMPARE_REALS(wvol, spacefill16::wvol5, spacefill16::eps2);
    COMPARE_REALS(wmean, spacefill16::wmean5, spacefill16::eps2);
    COMPARE_REALS(wgauss, spacefill16::wgauss5, spacefill16::eps2);
    final();
}

TEST_CASE("spacefill-17", "[spacefill][AMOEBA][pchloride]") {
    int argc = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/pchloride3.key",
        "../../test/testFiles/spacefill/pchloride3.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/pchloride4.key",
        "../../test/testFiles/spacefill/pchloride4.xyz",
        "1","Y","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/pchloride5.key",
        "../../test/testFiles/spacefill/pchloride5.xyz",
        "1","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);

    test = true;

    spacefill(argc, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill17::wsurf1, spacefill17::eps);
    COMPARE_REALS(wvol, spacefill17::wvol1, spacefill17::eps);
    COMPARE_REALS(wmean, spacefill17::wmean1, spacefill17::eps);
    COMPARE_REALS(wgauss, spacefill17::wgauss1, spacefill17::eps);
    final();

    spacefill(argc, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill17::wsurf2, spacefill17::eps2);
    COMPARE_REALS(wvol, spacefill17::wvol2, spacefill17::eps2);
    COMPARE_REALS(wmean, spacefill17::wmean2, spacefill17::eps2);
    COMPARE_REALS(wgauss, spacefill17::wgauss2, spacefill17::eps2);
    final();

    spacefill(argc, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill17::wsurf3, spacefill17::eps2);
    COMPARE_REALS(wvol, spacefill17::wvol3, spacefill17::eps2);
    COMPARE_REALS(wmean, spacefill17::wmean3, spacefill17::eps2);
    COMPARE_REALS(wgauss, spacefill17::wgauss3, spacefill17::eps2);
    final();
}

TEST_CASE("spacefill-18", "[spacefill][AMOEBA][crystal8k]") {
    // symmetric object; optional
    int argc = 9;
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/crystal8k_1.key",
        "../../test/testFiles/spacefill/crystal8k_1.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/crystal8k_2.key",
        "../../test/testFiles/spacefill/crystal8k_2.xyz",
        "1","Y","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/crystal8k_3.key",
        "../../test/testFiles/spacefill/crystal8k_3.xyz",
        "1","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    const char* strings4[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/crystal8k_4.key",
        "../../test/testFiles/spacefill/crystal8k_4.xyz",
        "1","Y","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill18::wsurf1, spacefill18::eps);
    COMPARE_REALS(wvol, spacefill18::wvol1, spacefill18::eps);
    COMPARE_REALS(wmean, spacefill18::wmean1, spacefill18::eps);
    COMPARE_REALS(wgauss, spacefill18::wgauss1, spacefill18::eps);
    final();

    spacefill(argc, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill18::wsurf2, spacefill18::eps);
    COMPARE_REALS(wvol, spacefill18::wvol2, spacefill18::eps);
    COMPARE_REALS(wmean, spacefill18::wmean2, spacefill18::eps);
    COMPARE_REALS(wgauss, spacefill18::wgauss2, spacefill18::eps);
    final();

    spacefill(argc, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill18::wsurf3, spacefill18::eps);
    COMPARE_REALS(wvol, spacefill18::wvol3, spacefill18::eps);
    COMPARE_REALS(wmean, spacefill18::wmean3, spacefill18::eps);
    COMPARE_REALS(wgauss, spacefill18::wgauss3, spacefill18::eps);
    final();

    spacefill(argc, argv4);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill18::wsurf4, spacefill18::eps);
    COMPARE_REALS(wvol, spacefill18::wvol4, spacefill18::eps);
    COMPARE_REALS(wmean, spacefill18::wmean4, spacefill18::eps);
    COMPARE_REALS(wgauss, spacefill18::wgauss4, spacefill18::eps);
    final();
}

TEST_CASE("spacefill-19", "[spacefill][AMOEBA][symmetry]") {
    int argc = 10;
    int n;

    // ermer: mirror symmetry
    const char* strings1[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/ermer.key",
        "../../test/testFiles/spacefill/ermer.xyz",
        "2","4.4","Y","N","Y","N"
    };
    char** argv1 = const_cast<char**>(strings1);

    // g9: mirror symmetry
    const char* strings2[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/g9.key",
        "../../test/testFiles/spacefill/g9.xyz",
        "2","1.4","Y","N","Y","N"
    };
    char** argv2 = const_cast<char**>(strings2);

    // cb8: center symmetry
    const char* strings3[] = {
        "spacefill",
        "-k", "../../test/testFiles/spacefill/cb8.key",
        "../../test/testFiles/spacefill/cb8.xyz",
        "2","1.4","Y","N","Y","N"
    };
    char** argv3 = const_cast<char**>(strings3);

    test = true;

    spacefill(argc, argv1);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill19::wsurf1, spacefill19::eps);
    COMPARE_REALS(wvol, spacefill19::wvol1, spacefill19::eps);
    n = 30;
    for (int i = 0; i < 3 * n; i++) {
        COMPARE_REALS(dsurf[i], spacefill19::dsurf1[i], spacefill19::eps);
        COMPARE_REALS(dvol[i], spacefill19::dvol1[i], spacefill19::eps);
    }
    final();

    spacefill(argc, argv2);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill19::wsurf2, spacefill19::eps);
    COMPARE_REALS(wvol, spacefill19::wvol2, spacefill19::eps);
    n = 26;
    for (int i = 0; i < 3 * n; i++) {
        COMPARE_REALS(dsurf[i], spacefill19::dsurf2[i], spacefill19::eps);
        COMPARE_REALS(dvol[i], spacefill19::dvol2[i], spacefill19::eps);
    }
    final();

    spacefill(argc, argv3);
    REQUIRE(alfmeth == AlfMethod::AlphaMol);
    COMPARE_REALS(wsurf, spacefill19::wsurf3, spacefill19::eps);
    COMPARE_REALS(wvol, spacefill19::wvol3, spacefill19::eps);
    n = 144;
    for (int i = 0; i < 3 * n; i++) {
        COMPARE_REALS(dsurf[i], spacefill19::dsurf3[i], spacefill19::eps);
        COMPARE_REALS(dvol[i], spacefill19::dvol3[i], spacefill19::eps);
    }
    final();
}
}
