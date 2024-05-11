cmake_minimum_required(VERSION 3.16.0)

set(POLMDQC_TEST
active.cpp
amoeba.cpp
amoebaMisc.cpp
amoebap.cpp
bounds.cpp
boyscheb.cpp
chkpole.cpp
chksymm.cpp
darray.cpp
getline.cpp
getnumb.cpp
getstring.cpp
gettext.cpp
getword.cpp
hartree.cpp
hippo.cpp
inertia.cpp
justify.cpp
lowcase.cpp
main.cpp
numeral.cpp
rotpole.cpp
sort.cpp
spacefill.cpp
spacefill2.cpp
spacefillT.cpp
torphase.cpp
trimhead.cpp
trimtext.cpp
upcase.cpp
)

add_executable(alltests ${POLMDQC_TEST})
target_link_libraries(alltests polmdqc openblas int2)
add_custom_target(test COMMAND alltests)
