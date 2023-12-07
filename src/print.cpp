////////////////////////////////////////////////
//                                            //
//  print.cpp  --  various printing routines  //
//                                            //
////////////////////////////////////////////////


#include "print.h"
#include <iostream>

namespace print
{
void printMatrix(const std::vector<std::vector<real>>& matrix, std::string header)
{
    int rowN = matrix.size();
    int colN = matrix[0].size();
    std::cout << header << std::endl;
    for (int i = 0; i < rowN; ++i)
    {
        printf("[");
        for (int j = 0; j < colN; ++j)
        {
            printf("%22.18f,", matrix[i][j]);
        }
        printf("],\n");
    }
}

void printVMatrix(const std::vector<real>& vmatrix, int rowN, int colN, std::string header)
{
    std::cout << header << std::endl;
    int N2 = rowN * colN;
    printf("[");
    for (int i = 0; i < N2; ++i)
    {
        printf("%22.18f,", vmatrix[i]);
        if (i == N2-1)
        {
            printf("],\n");
        }
        else if ((i+1) % rowN == 0)
            printf("],\n[");
    }
}
}
