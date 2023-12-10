// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "catch.hpp"
#include <vector>

#define COMPARE_REALS(val1, val2, eps) REQUIRE(fabs(val1 - val2) < eps)

#define COMPARE_VECTOR(vec1, vec2, eps)                  \
{                                                        \
    int n1 = vec1.size();                                \
    int n2 = vec2.size();                                \
    REQUIRE(n1 == n2);                                   \
    for (int i = 0; i < n1; i++) {                       \
        REQUIRE(fabs(vec1[i] - vec2[i]) < eps); \
    }                                                    \
}

#define COMPARE_MATRIX(mat1, mat2, eps)                            \
{                                                                  \
    int n1r = mat1.size();                                         \
    int n1c = (n1r > 0) ? mat1[0].size() : 0;                      \
    int n2r = mat2.size();                                         \
    int n2c = (n2r > 0) ? mat2[0].size() : 0;                      \
    REQUIRE(n1r == n2r);                                           \
    REQUIRE(n1c == n2c);                                           \
    for (int i = 0; i < n1r; i++) {                                \
        for (int j = 0; j < n1c; j++) {                            \
            REQUIRE(fabs(mat1[i][j] - mat2[i][j]) < eps); \
        }                                                          \
    }                                                              \
}
