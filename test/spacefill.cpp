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
    int argc = 7;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/water09.xyz",
        "1","Y","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    COMPARE_REALS(tsurf, spacefill1::tsurf, spacefill1::eps);
    COMPARE_REALS(tvol, spacefill1::tvol, spacefill1::eps);
    COMPARE_REALS(tmean, spacefill1::tmean, spacefill1::eps);
    COMPARE_REALS(tgauss, spacefill1::tgauss, spacefill1::eps);
    
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
    int argc = 7;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/water09.xyz",
        "1","N","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    COMPARE_REALS(tsurf, spacefill2::tsurf, spacefill2::eps);
    COMPARE_REALS(tvol, spacefill2::tvol, spacefill2::eps);
    COMPARE_REALS(tmean, spacefill2::tmean, spacefill2::eps);
    COMPARE_REALS(tgauss, spacefill2::tgauss, spacefill2::eps);

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
    int argc = 8;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","Y","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    COMPARE_REALS(tsurf, spacefill3::tsurf, spacefill3::eps);
    COMPARE_REALS(tvol, spacefill3::tvol, spacefill3::eps);
    COMPARE_REALS(tmean, spacefill3::tmean, spacefill3::eps);
    COMPARE_REALS(tgauss, spacefill3::tgauss, spacefill3::eps);

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
    int argc = 8;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","N","N","Y","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

    COMPARE_REALS(tsurf, spacefill4::tsurf, spacefill4::eps);
    COMPARE_REALS(tvol, spacefill4::tvol, spacefill4::eps);
    COMPARE_REALS(tmean, spacefill4::tmean, spacefill4::eps);
    COMPARE_REALS(tgauss, spacefill4::tgauss, spacefill4::eps);

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
    int argc = 9;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/water09.xyz",
        "2","1.4","Y","N","N","Y","1e-5"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);

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
    int argc1 = 7;
    const char* strings1[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 7;
    const char* strings2[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 8;
    const char* strings3[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 8;
    const char* strings4[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    COMPARE_REALS(tsurf, spacefill6::tsurf1, spacefill6::eps);
    COMPARE_REALS(tvol, spacefill6::tvol1, spacefill6::eps);
    COMPARE_REALS(tmean, spacefill6::tmean1, spacefill6::eps);
    COMPARE_REALS(tgauss, spacefill6::tgauss1, spacefill6::eps);
    final();

    spacefill(argc2, argv2);
    COMPARE_REALS(tsurf, spacefill6::tsurf2, spacefill6::eps);
    COMPARE_REALS(tvol, spacefill6::tvol2, spacefill6::eps);
    COMPARE_REALS(tmean, spacefill6::tmean2, spacefill6::eps);
    COMPARE_REALS(tgauss, spacefill6::tgauss2, spacefill6::eps);
    final();

    spacefill(argc3, argv3);
    COMPARE_REALS(tsurf, spacefill6::tsurf3, spacefill6::eps);
    COMPARE_REALS(tvol, spacefill6::tvol3, spacefill6::eps);
    COMPARE_REALS(tmean, spacefill6::tmean3, spacefill6::eps);
    COMPARE_REALS(tgauss, spacefill6::tgauss3, spacefill6::eps);
    final();

    spacefill(argc4, argv4);
    COMPARE_REALS(tsurf, spacefill6::tsurf4, spacefill6::eps);
    COMPARE_REALS(tvol, spacefill6::tvol4, spacefill6::eps);
    COMPARE_REALS(tmean, spacefill6::tmean4, spacefill6::eps);
    COMPARE_REALS(tgauss, spacefill6::tgauss4, spacefill6::eps);
    final();
}

TEST_CASE("spacefill-7", "[spacefill][AMOEBA][alatet]") {
    // inactive
    // 1
    int argc1 = 9;
    const char* strings1[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 9;
    const char* strings2[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 10;
    const char* strings3[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 10;
    const char* strings4[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/inactive.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    COMPARE_REALS(tsurf, spacefill7::tsurf1, spacefill7::eps);
    COMPARE_REALS(tvol, spacefill7::tvol1, spacefill7::eps);
    COMPARE_REALS(tmean, spacefill7::tmean1, spacefill7::eps);
    COMPARE_REALS(tgauss, spacefill7::tgauss1, spacefill7::eps);
    final();

    spacefill(argc2, argv2);
    COMPARE_REALS(tsurf, spacefill7::tsurf2, spacefill7::eps);
    COMPARE_REALS(tvol, spacefill7::tvol2, spacefill7::eps);
    COMPARE_REALS(tmean, spacefill7::tmean2, spacefill7::eps);
    COMPARE_REALS(tgauss, spacefill7::tgauss2, spacefill7::eps);
    final();

    spacefill(argc3, argv3);
    COMPARE_REALS(tsurf, spacefill7::tsurf3, spacefill7::eps);
    COMPARE_REALS(tvol, spacefill7::tvol3, spacefill7::eps);
    COMPARE_REALS(tmean, spacefill7::tmean3, spacefill7::eps);
    COMPARE_REALS(tgauss, spacefill7::tgauss3, spacefill7::eps);
    final();

    spacefill(argc4, argv4);
    COMPARE_REALS(tsurf, spacefill7::tsurf4, spacefill7::eps);
    COMPARE_REALS(tvol, spacefill7::tvol4, spacefill7::eps);
    COMPARE_REALS(tmean, spacefill7::tmean4, spacefill7::eps);
    COMPARE_REALS(tgauss, spacefill7::tgauss4, spacefill7::eps);
    final();
}

TEST_CASE("spacefill-8", "[spacefill][AMOEBA][alatet]") {
    // active
    // 1
    int argc1 = 9;
    const char* strings1[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv1 = const_cast<char**>(strings1);
    // 2
    int argc2 = 9;
    const char* strings2[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","N","N","N","N"
    };
    char** argv2 = const_cast<char**>(strings2);
    // 3
    int argc3 = 10;
    const char* strings3[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","Y","N","N","N"
    };
    char** argv3 = const_cast<char**>(strings3);
    // 4
    int argc4 = 10;
    const char* strings4[] = {
        "analyze",
        "-k", "../../test/testFiles/spacefill/active.key",
        "../../test/testFiles/spacefill/alatet.xyz",
        "2","1.4","N","N","N","N"
    };
    char** argv4 = const_cast<char**>(strings4);

    test = true;

    spacefill(argc1, argv1);
    COMPARE_REALS(tsurf, spacefill8::tsurf1, spacefill8::eps);
    COMPARE_REALS(tvol, spacefill8::tvol1, spacefill8::eps);
    COMPARE_REALS(tmean, spacefill8::tmean1, spacefill8::eps);
    COMPARE_REALS(tgauss, spacefill8::tgauss1, spacefill8::eps);
    final();

    spacefill(argc2, argv2);
    COMPARE_REALS(tsurf, spacefill8::tsurf2, spacefill8::eps);
    COMPARE_REALS(tvol, spacefill8::tvol2, spacefill8::eps);
    COMPARE_REALS(tmean, spacefill8::tmean2, spacefill8::eps);
    COMPARE_REALS(tgauss, spacefill8::tgauss2, spacefill8::eps);
    final();

    spacefill(argc3, argv3);
    COMPARE_REALS(tsurf, spacefill8::tsurf3, spacefill8::eps);
    COMPARE_REALS(tvol, spacefill8::tvol3, spacefill8::eps);
    COMPARE_REALS(tmean, spacefill8::tmean3, spacefill8::eps);
    COMPARE_REALS(tgauss, spacefill8::tgauss3, spacefill8::eps);
    final();

    spacefill(argc4, argv4);
    COMPARE_REALS(tsurf, spacefill8::tsurf4, spacefill8::eps);
    COMPARE_REALS(tvol, spacefill8::tvol4, spacefill8::eps);
    COMPARE_REALS(tmean, spacefill8::tmean4, spacefill8::eps);
    COMPARE_REALS(tgauss, spacefill8::tgauss4, spacefill8::eps);
    final();
}

TEST_CASE("spacefill-9", "[spacefill][AMOEBA][alatet]") {
    // atomic surface area, volume, mean curvature, and gaussian curvature
    int argc = 7;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","N","N"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);
    
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
    int argc = 8;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/spacefill/alatet.xyz",
        "1","Y","N","Y","Y","1e-5"
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    spacefill(argc, argv);
    
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
}
