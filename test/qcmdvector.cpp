#include "catch.hpp"
#include "qcmdvector.h"

TEST_CASE("qcmdvector") {
    Vector1D<int> qcmdv1;
    REQUIRE(qcmdv1.size() == 0);

    qcmdv1.allocate(10);
    REQUIRE(qcmdv1.size() == 10);

    qcmdv1.allocate(5);
    REQUIRE(qcmdv1.size() == 5);

    qcmdv1[0] = 0;
    qcmdv1[1] = 1;
    qcmdv1[2] = 2;
    qcmdv1[3] = 3;
    qcmdv1[4] = 4;
    int qcmdv10 = qcmdv1[0];
    int qcmdv11 = qcmdv1[1];
    int qcmdv12 = qcmdv1[2];
    int qcmdv13 = qcmdv1[3];
    int qcmdv14 = qcmdv1[4];
    REQUIRE(qcmdv1[0] == 0);
    REQUIRE(qcmdv1[1] == 1);
    REQUIRE(qcmdv1[2] == 2);
    REQUIRE(qcmdv1[3] == 3);
    REQUIRE(qcmdv1[4] == 4);
    REQUIRE(qcmdv10 == 0);
    REQUIRE(qcmdv11 == 1);
    REQUIRE(qcmdv12 == 2);
    REQUIRE(qcmdv13 == 3);
    REQUIRE(qcmdv14 == 4);

    Vector1D<int> qcmdv2(5);
    qcmdv2[0] = 4;
    qcmdv2[1] = 3;
    qcmdv2[2] = 2;
    qcmdv2[3] = 1;
    qcmdv2[4] = 0;
    REQUIRE(qcmdv2.size() == 5);
    REQUIRE(qcmdv2[0] == 4);
    REQUIRE(qcmdv2[1] == 3);
    REQUIRE(qcmdv2[2] == 2);
    REQUIRE(qcmdv2[3] == 1);
    REQUIRE(qcmdv2[4] == 0);
}
