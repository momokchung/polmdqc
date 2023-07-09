cmake_minimum_required(VERSION 3.16.0)

set(POLQCMD_TEST
gettext.cpp
main.cpp
upcase.cpp
)

add_executable(alltests ${POLQCMD_TEST})
target_link_libraries(alltests polqcmd)
add_custom_target(test COMMAND alltests)
