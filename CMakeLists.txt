cmake_minimum_required(VERSION 3.16.0)
project(tinkerqm)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the directory containing OpenBLAS libraries
link_directories(/Users/moseschung/OpenBLAS)

# Add the source files
set(SOURCES
source/tinkerqm.cpp
source/mechanic.cpp
source/katoms.cpp
source/kbasis.cpp
source/kgbs.cpp
source/ioUtils.cpp
source/stringUtils.cpp
source/mathUtils.cpp
source/config.cpp
source/energy.cpp
source/hartree.cpp
source/overlap.cpp
source/kprim.cpp
source/kinetic.cpp
source/boys.cpp
source/nuclear.cpp
source/init.cpp
source/nuclearRepulsion.cpp
source/kworker.cpp
)

# Add the executable
add_executable(tinkerqm ${SOURCES})

# Link the OpenBLAS library
target_link_libraries(tinkerqm openblas)
