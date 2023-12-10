#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "mpole.h"
#include "readjson.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("rotpole-1", "[AMOEBA][axetyp]") {
    // Tests Bisector, None, Z-Only, and Z-then-X
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/axetyp.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -4.0617143343279807;
    double emTest = -4.0617143343279807;
    double nemTest = 15;
    double eps = 1e-13;
    double epsR = 1e-16;

    std::string axetypJson = "../../test/testFiles/rotpole/axetyp.json";

    std::vector<std::vector<double>> refrpole = readMatrixFromJson(axetypJson);

    analyze(argc, argv);

    COMPARE_MATRIX(refrpole, rpole, epsR);
    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}

TEST_CASE("rotpole-2", "[AMOEBA][lysine_zbisect]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/lysine_zbisect.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double emTest = -68.364751107958028;
    double nemTest = 207;
    double eps = 1e-13;
    double epsR = 1e-16;

    std::string filename = "../../test/testFiles/rotpole/lysine_zbisect.json";
    std::vector<double> refrpole = readVectorFromJson(filename);

    analyze(argc, argv);

    COMPARE_VECTOR(refrpole, rpole[0], epsR);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}

TEST_CASE("rotpole-3", "[AMOEBA][lysine_3fold]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/lysine_3fold.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double emTest = -68.399778638605881;
    double nemTest = 207;
    double eps = 1e-13;
    double epsR = 1e-16;

    std::string filename = "../../test/testFiles/rotpole/lysine_3fold.json";
    std::vector<double> refrpole = readVectorFromJson(filename);

    analyze(argc, argv);

    COMPARE_VECTOR(refrpole, rpole[0], epsR);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}
}
