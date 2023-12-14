/*********************************************************************************
 *      A set of tools for sorting
 *********************************************************************************/

#pragma once
#include "vertex.h"
#include <vector>

/*********************************************************************************
  Class 
 *********************************************************************************/

class SortingTools
{
public:
    void valsort2(int a, int b, int *ia, int *ib, int *iswap);
    void valsort3(int a, int b, int c, int *ia, int *ib, int *ic, int *iswap);
    void valsort4(int a, int b, int c, int d, int *ia, int *ib, int *ic, int *id, int *iswap);
    void valsort5(int a, int b, int c, int d, int e, 
        int *ia, int *ib, int *ic, int *id, int *ie, int *iswap);
    void isort_indx(int *list, int *idx, int *nswap, int n);
    void isort_swap(int *list, int *nswap, int n);
    void isort4_swap(int *a, int *b, int *c, int *d, int *nswap);
    void sort4_sign(int *list, int *idx, int *nswap, int n);
};
