cmake_minimum_required(VERSION 3.16.0)
project(polqcmd)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add the directory containing OpenBLAS libraries
link_directories(/Users/moseschung/OpenBLAS)

# Add the directories containing header files
include_directories(include)

# Add the source files
include (cmakesrc.txt)
add_library (polqcmd STATIC "${POLQCMD_SRC}")

# Add the tinkerqm executable
add_executable(tinkerqm src/tinkerqm.cpp)

# Link the tinkerqm libraries
target_link_libraries(tinkerqm polqcmd openblas)

# Add the analyze executable
add_executable(analyze src/analyze.cpp)

# Link the polqcmd library
target_link_libraries(analyze polqcmd)

# Include tests directory
add_subdirectory(test)
