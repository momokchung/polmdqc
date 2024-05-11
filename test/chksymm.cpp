#include "atomid.h"
#include "atoms.h"
#include "chksymm.h"
#include "files.h"
#include "final.h"
#include "getcart.h"
#include "inform.h"
#include "initial.h"
#include "mechanic.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("chksymm-1", "[chksymm][None]") {
    int argc = 2;
    const char* strings[] = {
        "dummy",
        "../../test/testFiles/chksymm/alatet.xyz",
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    // initialize
    initial(argc, argv);
    getcart(ffile);
    mechanic();

    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);

    REQUIRE(symtyp == SymTyp::None);

    ffile.close();
    final();
}

TEST_CASE("chksymm-2", "[chksymm][Single]") {
    int argc = 2;
    const char* strings[] = {
        "dummy",
        "../../test/testFiles/chksymm/chloride.xyz",
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    // initialize
    initial(argc, argv);
    getcart(ffile);
    mechanic();

    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);

    REQUIRE(symtyp == SymTyp::Single);

    ffile.close();
    final();
}

TEST_CASE("chksymm-3", "[chksymm][Linear]") {
    int argc1 = 4;
    const char* strings1[] = {
        "dummy",
        "-k", "../../test/testFiles/chksymm/chloride.key",
        "../../test/testFiles/chksymm/lchloride2.xyz",
    };
    char** argv1 = const_cast<char**>(strings1);

    int argc2 = 4;
    const char* strings2[] = {
        "dummy",
        "-k", "../../test/testFiles/chksymm/chloride.key",
        "../../test/testFiles/chksymm/lchloride3.xyz",
    };
    char** argv2 = const_cast<char**>(strings2);

    test = true;

    initial(argc1, argv1);
    getcart(ffile);
    mechanic();
    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Linear);
    ffile.close();
    final();

    initial(argc2, argv2);
    getcart(ffile);
    mechanic();
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Linear);
    ffile.close();
    final();
}

TEST_CASE("chksymm-4", "[chksymm][Planar]") {
    int argc1 = 4;
    const char* strings1[] = {
        "dummy",
        "-k", "../../test/testFiles/chksymm/chloride.key",
        "../../test/testFiles/chksymm/pchloride3.xyz",
    };
    char** argv1 = const_cast<char**>(strings1);

    int argc2 = 4;
    const char* strings2[] = {
        "dummy",
        "-k", "../../test/testFiles/chksymm/chloride.key",
        "../../test/testFiles/chksymm/pchloride4.xyz",
    };
    char** argv2 = const_cast<char**>(strings2);

    test = true;

    initial(argc1, argv1);
    getcart(ffile);
    mechanic();
    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Planar);
    ffile.close();
    final();

    initial(argc2, argv2);
    getcart(ffile);
    mechanic();
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Planar);
    ffile.close();
    final();
}

TEST_CASE("chksymm-5", "[chksymm][Mirror]") {
    int argc = 2;
    const char* strings[] = {
        "dummy",
        "../../test/testFiles/chksymm/g9.xyz",
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    initial(argc, argv);
    getcart(ffile);
    mechanic();
    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Mirror);
    ffile.close();
    final();
}

TEST_CASE("chksymm-6", "[chksymm][Center]") {
    int argc = 2;
    const char* strings[] = {
        "dummy",
        "../../test/testFiles/chksymm/cb8.xyz",
    };
    char** argv = const_cast<char**>(strings);

    test = true;

    initial(argc, argv);
    getcart(ffile);
    mechanic();
    SymTyp symtyp;
    chksymm(n, mass.ptr(), x.ptr(), y.ptr(), z.ptr(), symtyp);
    REQUIRE(symtyp == SymTyp::Center);
    ffile.close();
    final();
}
}