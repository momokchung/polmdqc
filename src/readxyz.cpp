////////////////////////////////////////////////////////
//                                                    //
//  readxyz.cpp  --  input of XYZ-format coordinates  //
//                                                    //
////////////////////////////////////////////////////////

// "readxyz" gets a set of Cartesian coordinates from
// an external disk file


#include "atomid.h"
#include "atoms.h"
#include "boxes.h"
#include "couple.h"
#include "fatal.h"
#include "files.h"
#include "getline.h"
#include "gettext.h"
#include "getword.h"
#include "inform.h"
#include "inquire.h"
#include "lattice.h"
#include "readxyz.h"
#include "titles.h"
#include "trimtext.h"
#include "unitcell.h"
#include <sstream>

void readxyz(std::string& xyzfile)
{
    int nmax;
    int next,size;
    std::vector<int> list;
    bool exist;
    bool quit,reorder;
    bool clash;
    std::string record;
    std::string string;
    std::string line;
    std::istringstream iss;

    // initialize the total number of atoms in the system
    n = 0;

    // open the input file
    exist = inquire(xyzfile);
    if (!exist) {
        printf("\n READXYZ  --  Unable to Find the Cartesian Coordinates File\n");
        fatal();
    }
    std::ifstream file (xyzfile);

    // read first line and return if already at end of file
    quit = false;
    informAbort = true;
    size = 0;
    while (std::getline(file, record)) {
        iss.str(record);
        std::string firstWord;
        iss >> firstWord;
        size = firstWord.length();
        if (size != 0) {
            break;
        }
    }
    if (size == 0) {
        goto label_80;
    }
    informAbort = false;
    quit = true;

    // parse the title line to get the number of atoms
    next = 0;
    gettext(record, string, next);
    iss.str(string);
    if (!(iss >> n)) goto label_80;

    // extract the title and determine its length
    string = record.substr(next);
    line = getline(string);
    ltitle = line.length();
    title = line;
    if (ltitle == 0) title = " ";

    // check for too few or too many total atoms in the file
    if (n <= 0) {
        printf("\n READXYZ  --  The Coordinate File Does Not Contain Any Atoms\n");
        fatal();
    }
    else if (n > maxatm) {
        printf("\n  READXYZ  --  The Maximum of %d Atoms has been Exceeded\n", maxatm);
        fatal();
    }

    // initialize coordinates and connectivities for each atom
    for (int i = 0; i < n; i++) {
        tag[i] = 0;
        name[i] = "   ";
        x[i] = 0.;
        y[i] = 0.;
        z[i] = 0.;
        type[i] = 0;
        n12[i] = 0;
        for (int j = 0; j < maxval; j++) {
            i12[i][j] = 0;
        }
    }

    // read the coordinates and connectivities for each atom
    for (int i = 0; i < n; i++) {
        size = 0;
        while (size == 0) {
            unitcell();
            std::getline(file, record);
            size = trimtext(record);
            if (i == 0) {
                next = 0;
                getword (record, name[i], next);
                if (name[i].length() == 0) {
                    iss.str(record);
                    if ((iss >> xbox >> ybox >> zbox >> alpha >> beta >> gamma)) size = 0;
                }
            }
            lattice();
        }
//          read (record,*,err=80,end=80)  tag(i)
//          next = 1
//          call getword (record,name(i),next)
//          string = record(next:240)
//          read (string,*,err=70,end=70)  x(i),y(i),z(i),type(i),
//      &                                  (i12(j,i),j=1,maxval)
//    70    continue
    }
//       quit = .false.
//    80 continue
//       if (.not. opened)  close (unit=ixyz)


    label_80:

    printf("");
}

// c
// c     read the coordinates and connectivities for each atom
// c
//       do i = 1, n
//          size = 0
//          do while (size .eq. 0)
//             call unitcell
//             read (ixyz,50,err=80,end=80)  record
//    50       format (a240)
//             size = trimtext (record)
//             if (i .eq. 1) then
//                next = 1
//                call getword (record,name(i),next)
//                if (name(i) .ne. '   ')  goto 60
//                read (record,*,err=60,end=60)  xbox,ybox,zbox,
//      &                                        alpha,beta,gamma
//                size = 0
//             end if
//    60       continue
//             call lattice
//          end do
//          read (record,*,err=80,end=80)  tag(i)
//          next = 1
//          call getword (record,name(i),next)
//          string = record(next:240)
//          read (string,*,err=70,end=70)  x(i),y(i),z(i),type(i),
//      &                                  (i12(j,i),j=1,maxval)
//    70    continue
//       end do
//       quit = .false.
//    80 continue
//       if (.not. opened)  close (unit=ixyz)
// c
// c     an error occurred in reading the coordinate file
// c
//       if (quit) then
//          write (iout,90)  i
//    90    format (/,' READXYZ  --  Error in Coordinate File at Atom',i9)
//          call fatal
//       end if
// c
// c     for each atom, count and sort its attached atoms
// c
//       do i = 1, n
//          do j = maxval, 1, -1
//             if (i12(j,i) .ne. 0) then
//                n12(i) = j
//                goto 100
//             end if
//          end do
//   100    continue
//          call sort (n12(i),i12(1,i))
//       end do
// c
// c     perform dynamic allocation of some local arrays
// c
//       nmax = 0
//       do i = 1, n
//          nmax = max(tag(i),nmax)
//          do j = 1, n12(i)
//             nmax = max(i12(j,i),nmax)
//          end do
//       end do
//       allocate (list(nmax))
// c
// c     check for scrambled atom order and attempt to renumber
// c
//       reorder = .false.
//       do i = 1, n
//          list(tag(i)) = i
//          if (tag(i) .ne. i)  reorder = .true.
//       end do
//       if (reorder) then
//          write (iout,110)
//   110    format (/,' READXYZ  --  Atom Labels not Sequential,',
//      &              ' Attempting to Renumber')
//          do i = 1, n
//             tag(i) = i
//             do j = 1, n12(i)
//                i12(j,i) = list(i12(j,i))
//             end do
//             call sort (n12(i),i12(1,i))
//          end do
//       end if
// c
// c     perform deallocation of some local arrays
// c
//       deallocate (list)
// c
// c     check for atom pairs with identical coordinates
// c
//       clash = .false.
//       if (n .le. 10000)  call chkxyz (clash)
// c
// c     make sure all atom connectivities are bidirectional
// c
//       do i = 1, n
//          do j = 1, n12(i)
//             k = i12(j,i)
//             do m = 1, n12(k)
//                if (i12(m,k) .eq. i)  goto 130
//             end do
//             write (iout,120)  k,i
//   120       format (/,' READXYZ  --  Check Connection of Atoms',
//      &                 i9,' and',i9)
//             call fatal
//   130       continue
//          end do
//       end do
//       return
//       end
