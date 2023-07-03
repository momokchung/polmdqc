////////////////////////////////////
//                                //
//  kworker.h  --  worker arrays  //
//                                //
////////////////////////////////////


#include "init.h"

namespace worker
{
extern std::vector<real> workerN_1;
extern std::vector<real> workerN_2;
extern std::vector<real> workerN2_1;
extern std::vector<real> workerN2_2;

void allocateWorker();
}
