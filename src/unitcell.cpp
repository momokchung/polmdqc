// //////////////////////////////////////////////////////////
// //                                                      //
// //  unitcell.cpp  --  get periodic boundary conditions  //
// //                                                      //
// //////////////////////////////////////////////////////////

// // "unitcell" gets the periodic boundary box size and related
// // values from an external keyword file


// #include "bound.h"
// #include "boxes.h"
// #include "gettext.h"
// #include "keys.h"
// #include "math.h"
// #include "unitcell.h"
// #include "upcase.h"

// void unitcell()
// {
//     // set the default values for periodic boundary conditions
//     use_bounds = false;
//     use_replica = false;

//     // set the default values for the unit cell variables
//     orthogonal = false;
//     monoclinic = false;
//     triclinic = false;
//     octahedron = false;
//     dodecadron = false;
//     nonprism = false;
//     nosymm = false;
//     spacegrp = "          ";

//     // get keywords containing crystal lattice dimensions
//     for (int i = 0; i < nkey; i++) {
//         int next = 1;
//         std::string record = keyline[i];
//         std::string keyword;
//         gettext(record, keyword, next);
//         upcase(keyword);
//         std::string string = record.substr(next);
//         std::string nextString;
//         if (keyword == "X-AXIS") {
//             if (xbox == 0.) {
//                 gettext(string,nextString,next);
//                 try {
//                     xbox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "Y-AXIS") {
//             if (ybox == 0.)  {
//                 gettext(string,nextString,next);
//                 try {
//                     ybox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "Z-AXIS") {
//             if (zbox == 0.)  {
//                 gettext(string,nextString,next);
//                 try {
//                     zbox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "A-AXIS") {
//             if (xbox == 0.)  {
//                 gettext(string,nextString,next);
//                 try {
//                     xbox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "B-AXIS") {
//             if (ybox == 0.)  {
//                 gettext(string,nextString,next);
//                 try {
//                     ybox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "C-AXIS") {
//             if (zbox == 0.)  {
//                 gettext(string,nextString,next);
//                 try {
//                     zbox = std::stod(keyword);
//                 } catch (const std::exception& e) {}
//             }
//         }
//         else if (keyword == "ALPHA") {
//             if (alpha == 0.)  read (string,*,err=10,end=10)  alpha
//         }
//         else if (keyword == "BETA") {
//             if (beta == 0.)  read (string,*,err=10,end=10)  beta
//         }
//         else if (keyword == "GAMMA") {
//             if (gamma == 0.)  read (string,*,err=10,end=10)  gamma
//         }
//         else if (keyword == "OCTAHEDRON") {
//             octahedron = true;
//         }
//         else if (keyword == "DODECAHEDRON") {
//             dodecadron = true;
//         }
//         else if (keyword == "NOSYMMETRY") {
//             nosymm = true;
//         }
//         else if (keyword(1:11) == "SPACEGROUP") {
//             call getword (record,spacegrp,next)
//         }
//     }
// }

// c

// c
// c     use periodic boundary conditions if a cell was defined
// c
//       boxmax = max(xbox,ybox,zbox)
//       if (boxmax .ne. 0.0d0)  use_bounds = .true.
// c
// c     set unspecified periodic boundary box lengths and angles
// c
//       if (use_bounds) then
//          if (xbox .eq. 0.0d0)  xbox = boxmax
//          if (ybox .eq. 0.0d0)  ybox = boxmax
//          if (zbox .eq. 0.0d0)  zbox = boxmax
//          if (alpha .eq. 0.0d0)  alpha = 90.0d0
//          if (beta .eq. 0.0d0)  beta = 90.0d0
//          if (gamma .eq. 0.0d0)  gamma = 90.0d0
// c
// c     determine the general periodic boundary lattice type
// c
//          if (nosymm) then
//             triclinic = .true.
//          else if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
//      &               .and. gamma.eq.90.0d0) then
//             orthogonal = .true.
//          else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
//             monoclinic = .true.
//          else
//             triclinic = .true.
//          end if
//       end if
// c
// c     set lattice values for non-prism periodic boundaries
// c
//       if (octahedron .or. dodecadron) then
//          orthogonal = .false.
//          monoclinic = .false.
//          triclinic = .false.
//          nonprism = .true.
//          ybox = xbox
//          if (octahedron) then
//             zbox = xbox
//          else if (dodecadron) then
//             zbox = xbox * root2
//          end if         
//          alpha = 90.0d0
//          beta = 90.0d0
//          gamma = 90.0d0
//       end if
//       return
//       end
