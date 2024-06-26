cmake_minimum_required(VERSION 3.16.0)
set(CMAKE_C_COMPILER "/opt/homebrew/Cellar/gcc/13.2.0/bin/gcc-13")
set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/gcc/13.2.0/bin/g++-13")
project(polmdqc)
# set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_OSX_SYSROOT "")

set(CMAKE_C_COMPILER "/opt/homebrew/Cellar/gcc/13.2.0/bin/gcc-13")
set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/gcc/13.2.0/bin/g++-13")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-functions -march=native -ffast-math -fopenmp")
set(CMAKE_CXX_FLAGS "-O3 -finline-functions -ffast-math -fopenmp")
# Debug
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wunused-parameter -Wunused-variable -Wunused-function -fopenmp")

# Link directories
# libint2
link_directories(/Users/moseschung/libint-2.7.2)
# OpenBLAS
link_directories(/Users/moseschung/OpenBLAS)
# gmp
link_directories(/opt/homebrew/Cellar/gmp/6.3.0/lib)

# Include directories
include_directories(/usr/local/include)
include_directories(include)
include_directories(include/md)
include_directories(include/oqc)
include_directories(include/qc)
# boost
include_directories(/opt/homebrew/Cellar/boost/1.83.0/include)
# gmp
include_directories(/opt/homebrew/Cellar/gmp/6.3.0/include)

# Include
include(cmakesrc.txt)

# Add library
add_library (polmdqc STATIC "${POLMDQC_SRC}")

# Add executable
# analyze
add_executable(analyze src/md/analyzeExe.cpp)
# testgrad
add_executable(testgrad src/md/testgradExe.cpp)
# spacefill
add_executable(spacefill src/md/spacefillExe.cpp)
# tinkerqm
add_executable(tinkerqm src/oqc/tinkerqmExe.cpp)
# analyzeqm
add_executable(analyzeqm src/qc/analyzeqmExe.cpp)
# valeev
add_executable(valeev1 src/oqc/valeev1.cpp)

# Target link libraries
# analyze
target_link_libraries(analyze polmdqc)
# testgrad
target_link_libraries(testgrad polmdqc)
# spacefill
target_link_libraries(spacefill polmdqc gmp)
# tinkerqm
target_link_libraries(tinkerqm polmdqc openblas int2)
# analyzeqm
target_link_libraries(analyzeqm polmdqc openblas)
# valeev
target_link_libraries(valeev1 polmdqc int2)

# Add subdirectory
add_subdirectory(test)
