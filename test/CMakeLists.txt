cmake_minimum_required(VERSION 3.16.0)

set(POLQCMD_TEST
analyze.cpp
getline.cpp
getnumb.cpp
getstring.cpp
gettext.cpp
getword.cpp
hartree.cpp
justify.cpp
lowcase.cpp
main.cpp
numeral.cpp
qcmdvector.cpp
rotpole.cpp
sort.cpp
torphase.cpp
trimhead.cpp
trimtext.cpp
upcase.cpp
)

add_executable(alltests ${POLQCMD_TEST})
target_link_libraries(alltests polqcmd openblas int2)
add_custom_target(test COMMAND alltests)
