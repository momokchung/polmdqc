//////////////////////////////////////
//                                  //
//  kworker.cpp  --  worker arrays  //
//                                  //
//////////////////////////////////////


#include "kbasis.h"
#include "kworker.h"
#include <vector>

namespace worker
{
// workerN_1     worker1 of basis length
// workerN_2     worker2 of basis length
// workerN2_1    worker1 of basis length squared
// workerN2_2    worker2 of basis length squared

std::vector<real> workerN_1;
std::vector<real> workerN_2;
std::vector<real> workerN2_1;
std::vector<real> workerN2_2;


///////////////////////////////////////////////////////
//                                                   //
//  void allocateWorker  --  allocate worker arrays  //
//                                                   //
///////////////////////////////////////////////////////

void allocateWorker()
{
    int N = basis::N;
    int N2 = N * N;

    workerN_1.resize(0);
    workerN_2.resize(0);
    workerN2_1.resize(0);
    workerN2_2.resize(0);

    workerN_1.reserve(N);
    workerN_2.reserve(N);
    workerN2_1.reserve(N2);
    workerN2_2.reserve(N2);
}
}
