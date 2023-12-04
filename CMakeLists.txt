cmake_minimum_required(VERSION 3.16.0)
project(polmdqc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-functions -march=native -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-functions -ffast-math")

# Add the directory containing OpenBLAS libraries
link_directories(/Users/moseschung/OpenBLAS)

# Add libint2 package
link_directories(/Users/moseschung/libint-2.7.2)
find_package(Libint2)
include_directories(/usr/local/include)

# Include boost directory
include_directories(/opt/homebrew/Cellar/boost/1.83.0/include)

# Add the directories containing header files
include_directories(include)

# Add the source files
include (cmakesrc.txt)
add_library (polmdqc STATIC "${POLMDQC_SRC}")

# Add the tinkerqm executable
add_executable(tinkerqm src/tinkerqmExe.cpp)

# Link the tinkerqm libraries
target_link_libraries(tinkerqm polmdqc openblas int2)

# Add the analyze executable
add_executable(analyze src/analyzeExe.cpp)

# Link the polmdqc library
target_link_libraries(analyze polmdqc)

# Add the alphamol executable
add_executable(alphamol src/alphamolExe.cpp)

# Link the polmdqc library
target_link_libraries(alphamol polmdqc)

# Include tests directory
add_subdirectory(test)

# Link the tinkerqm libraries
add_executable(valeev1 src/valeev1.cpp)
target_link_libraries(valeev1 polmdqc int2)
