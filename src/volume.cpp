///////////////////////////////////////////////////////
//                                                   //
//  volume.cpp  --  compute surface area and volume  //
//                                                   //
///////////////////////////////////////////////////////

// "volume" computes surface area and volume using GaussVol procedure


#include "active.h"
#include "atomid.h"
#include "atoms.h"
#include "field.h"
#include "files.h"
#include "final.h"
#include "gaussvol.h"
#include "gettext.h"
#include "getxyz.h"
#include "inform.h"
#include "initial.h"
#include "katom.h"
#include "kvdw.h"
#include "kvdws.h"
#include "mathConst.h"
#include "nextarg.h"
#include "ptable.h"
#include "readxyz.h"
#include "suffix.h"
#include "upcase.h"
#include "usage.h"
#include "volume.h"
#include <cmath>
#include <iostream>
#include <vector>
      
void volume(int argc, char** argv)
{
    int next;
    int frame;
    std::vector<int> ishydrogen;
    double gv_volume,gv_area,gv_energy,gv_volume2;
    double probe,exclude,rad_offset;
    std::vector<double> posx;
    std::vector<double> posy;
    std::vector<double> posz;
    std::vector<double> radius;
    std::vector<double> vol;
    std::vector<double> radius2;
    std::vector<double> vol2;
    std::vector<double> drx;
    std::vector<double> dry;
    std::vector<double> drz;
    std::vector<double> dv;
    std::vector<double> free_volume;
    std::vector<double> self_volume;
    bool exist,query,use_hydrogen;
    std::string answer;
    std::string xyzfile;
    std::string record;
    std::string string;

    // get the Cartesian coordinates for the system
    initial(argc, argv);
    getxyz();

    // determine the atoms to be used in computation;
    // radii can be changed via the keyword mechanism
    field();
    active();
    katom();
    kvdw();

    // initialize random numbers and turn on extra printing
    verbose = false;
    debug = true;

    // decide whether to include hydrogens in the calculation
    nextarg(answer,exist);
    if (!exist) {
        printf("\n Include the Hydrogen Atoms in Computation [N] :  ");
        std::getline(std::cin, record);
        next = 0;
        gettext(record,answer,next);
    }
    use_hydrogen = true;
    upcase(answer);
    if (answer != "Y") {
        use_hydrogen = false;
        for (int i = 0; i < n; i++) {
            if (atomic[i] == 1) use[i] = false;
        }
    }

    // decide whether to provide full output for large systems
    if (n > 100) {
        debug = false;
        nextarg(answer,exist);
        if (!exist) {
            printf("\n Output the Surface Area of Individual Atoms [N] :  ");
            std::getline(std::cin, record);
            next = 0;
            gettext(record,answer,next);
        }
        upcase(answer);
        if (answer == "Y") debug = true;
    }

    // perform dynamic allocation of some local arrays
    radius.resize(n);
    radius2.resize(n);
    ishydrogen.resize(n, 0);
    posx.resize(n);
    posy.resize(n);
    posz.resize(n);
    vol.resize(n);
    vol2.resize(n);
    drx.resize(n);
    dry.resize(n);
    drz.resize(n);
    dv.resize(n);
    free_volume.resize(n);
    self_volume.resize(n);
    for (int i = 0; i < n; i++) {
        if (atomic[i] == 1 and !use_hydrogen) ishydrogen[i] = 1;
        posx[i] = x[i];
        posy[i] = y[i];
        posz[i] = z[i];
    }

    // set atomic radii based on force field or Bondi values
    rad_offset = 0.005;
    // rad_offset = 0.001;
    for (int i = 0; i < n; i++) {
        radius[i] = rad[atomClass[i]];
        radius2[i] = rad[atomClass[i]] + rad_offset;
        // radius[i] = rad[atomClass[i]] / twosix;
        // radius[i] = vdwrad(atomic[i]);
    }
    constexpr double four_thirds_pi = 4./3.*pi;
    for (int i = 0; i < n; i++) {
        if (use[i]) {
            vol[i] = four_thirds_pi*std::pow(radius[i],3);
            vol2[i] = four_thirds_pi*std::pow(radius2[i],3);
        }
        else {
            vol[i] = 0.;
            vol2[i] = 0.;
        }
    }

    // reopen the coordinates file and read the first structure
    frame = 0;
    xyzfile = filename;
    suffix(xyzfile,"xyz","old");
    std::ifstream ffile(xyzfile);
    readxyz(ffile);

    // GaussVol procedure
    GaussVol gvol(n, ishydrogen);
    gvol.setRadii(radius);
    gvol.setVolumes(vol);
    gvol.compute_tree(posx, posy, posz);
    // gvol.print_tree();
    gvol.compute_volume(posx, posy, posz, gv_volume, gv_energy, drx, dry, drz, dv, free_volume, self_volume);
    // write(*,*) "GaussVol Volume:  ", volume
    gvol.setRadii(radius2);
    gvol.setVolumes(vol2);
    gvol.rescan_tree_volumes(posx, posy, posz);
    gvol.compute_volume(posx, posy, posz, gv_volume2, gv_energy, drx, dry, drz, dv, free_volume, self_volume);
    gv_area = (gv_volume2 - gv_volume)/rad_offset;
    printf(" GaussVol Volume: %20.16e\n", gv_volume);
    printf(" GaussVol Volume2: %20.16e\n", gv_volume2);
    printf(" GaussVol Area: %20.16e\n", gv_area);

// c     get area and volume for successive coordinate structures
// c
//       do while (.not. abort)
//          frame = frame + 1
//          if (frame .gt. 1) then
//             write (iout,120)  frame
//   120       format (/,' Area and Volume for Archive Structure :',5x,i8)
//          end if
// c
// c     use the Connolly routines to find the area and volume
// c
//          call connolly (volume,area,radius,probe,exclude)
// c
// c     print out the values of the total surface area and volume
// c
//          if (mode .eq. 1) then
//             write (iout,130)
//   130       format (/,' Van der Waals Surface Area and Volume :')
//          else if (mode .eq. 2) then
//             write (iout,140)
//   140       format (/,' Accessible Surface Area and Excluded Volume :')
//          else if (mode .eq. 3) then
//             write (iout,150)
//   150       format (/,' Contact-Reentrant Surface Area and Volume :')
//          end if
//          write (iout,160)  area
//   160    format (/,' Total Area :',f20.3,' Square Angstroms')
//          write (iout,170)  volume
//   170    format (' Total Volume :',f18.3,' Cubic Angstroms')
// c
// c     attempt to read next structure from the coordinate file
// c
//          call readxyz (ixyz)
//       end do

    // perform any final tasks before program exit
    ffile.close();
    final();
}
