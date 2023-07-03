cmake_minimum_required(VERSION 3.16.0)
project(polqcmd)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the directory containing OpenBLAS libraries
link_directories(/Users/moseschung/OpenBLAS)

# Add the source files
add_library(polqcmd
source/boys.cpp
source/config.cpp
source/energy.cpp
source/hartree.cpp
source/init.cpp
source/initial.cpp
source/ioUtils.cpp
source/katoms.cpp
source/kbasis.cpp
source/kgbs.cpp
source/kinetic.cpp
source/kprim.cpp
source/kworker.cpp
source/mathUtils.cpp
source/mechanic.cpp
source/nuclear.cpp
source/nuclearRepulsion.cpp
source/overlap.cpp
source/promo.cpp
source/stringUtils.cpp
)

# Add the tinkerqm executable
add_executable(tinkerqm source/tinkerqm.cpp)

# Link the tinkerqm libraries
target_link_libraries(tinkerqm polqcmd openblas)

# Add the analyze executable
add_executable(analyze source/analyze.cpp)

# Link the polqcmd library
target_link_libraries(analyze polqcmd)

