cmake_minimum_required(VERSION 3.16.0)
project(polqcmd)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the directory containing OpenBLAS libraries
link_directories(/Users/moseschung/OpenBLAS)

# Add the directories containing header files
include_directories(include)

# Add the source files
add_library(polqcmd
src/boys.cpp
src/command.cpp
src/config.cpp
src/energy.cpp
src/hartree.cpp
src/init.cpp
src/initatom.cpp
src/initial.cpp
src/ioUtils.cpp
src/katoms.cpp
src/kbasis.cpp
src/kgbs.cpp
src/kinetic.cpp
src/kprim.cpp
src/kworker.cpp
src/mathUtils.cpp
src/mechanic.cpp
src/nuclear.cpp
src/nuclearRepulsion.cpp
src/overlap.cpp
src/promo.cpp
src/stringUtils.cpp
)

# Add the tinkerqm executable
add_executable(tinkerqm src/tinkerqm.cpp)

# Link the tinkerqm libraries
target_link_libraries(tinkerqm polqcmd openblas)

# Add the analyze executable
add_executable(analyze src/analyze.cpp)

# Link the polqcmd library
target_link_libraries(analyze polqcmd)

