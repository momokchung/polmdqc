#include "kinetic.h"
#include "nuclear.h"
#include "overlap.h"
#include "readjson.h"
#include "testrt.h"
#include "tinkerqm.h"

TEST_CASE("tinkerqm-1", "[HartreeFock][water_3-21g]") {
    int argc = 3;
    const char* strings[] = {
        "tinkerqm",
        "../../test/testFiles/water_3-21g/water_3-21g.xyz",
    };
    char** argv = const_cast<char**>(strings);

    tinkerqm(argc, argv);

    std::string overlapJson = "../../test/testFiles/water_3-21g/overlap.json";
    std::string kineticJson = "../../test/testFiles/water_3-21g/kinetic.json";
    std::string nuclearJson = "../../test/testFiles/water_3-21g/nuclear.json";

    std::vector<std::vector<double>> refcartS = readMatrixFromJson(overlapJson);
    std::vector<std::vector<double>> refcartKE = readMatrixFromJson(kineticJson);
    std::vector<std::vector<double>> refcartNE = readMatrixFromJson(nuclearJson);

    double epsS = 1e-15;
    double epsKE = 1e-14;
    double epsNE = 1e-13;
    COMPARE_MATRIX(refcartS, overlap::cartS, epsS);
    COMPARE_MATRIX(refcartKE, kinetic::cartKE, epsKE);
    COMPARE_MATRIX(refcartNE, nuclear::cartNE, epsNE);
}

TEST_CASE("tinkerqm-2", "[HartreeFock][water_aug-cc-pvtz]") {
    int argc = 3;
    const char* strings[] = {
        "tinkerqm",
        "../../test/testFiles/water_aug-cc-pvtz/water_aug-cc-pvtz.xyz",
    };
    char** argv = const_cast<char**>(strings);

    tinkerqm(argc, argv);

    std::string overlapJson = "../../test/testFiles/water_aug-cc-pvtz/overlap.json";
    std::string kineticJson = "../../test/testFiles/water_aug-cc-pvtz/kinetic.json";
    std::string nuclearJson = "../../test/testFiles/water_aug-cc-pvtz/nuclear.json";

    std::vector<std::vector<double>> refsphS = readMatrixFromJson(overlapJson);
    std::vector<std::vector<double>> refsphKE = readMatrixFromJson(kineticJson);
    std::vector<std::vector<double>> refsphNE = readMatrixFromJson(nuclearJson);

    double epsS = 1e-14;
    double epsKE = 1e-13;
    double epsNE = 1e-13;
    COMPARE_MATRIX(refsphS, overlap::sphS, epsS);
    COMPARE_MATRIX(refsphKE, kinetic::sphKE, epsKE);
    COMPARE_MATRIX(refsphNE, nuclear::sphNE, epsNE);
}
