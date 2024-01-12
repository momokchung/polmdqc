// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "bound.h"
#include "cutoffs.h"
#include "gettext.h"
#include "hescut.h"
#include "keys.h"
#include "neigh.h"
#include "polpot.h"
#include "mdqclimits.h"
#include "tarray.h"
#include "upcase.h"
#include <cmath>
#include <sstream>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  cutoffs  --  distance cutoffs & neighbor lists  //
//                                                  //
//////////////////////////////////////////////////////

// "cutoffs" initializes and stores spherical energy cutoff
// distance windows, Hessian element and Ewald sum cutoffs,
// and allocates pairwise neighbor lists

void cutoffs()
{
    int next;
    int limit;
    real big,value;
    bool truncate;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // set defaults for spherical energy cutoff distances
    big = 1.0e12;
    if (use_bounds) {
        vdwcut = 9.0;
        dispcut = 9.0;
        chgcut = 9.0;
        dplcut = 9.0;
        mpolecut = 9.0;
    }
    else {
        vdwcut = big;
        dispcut = big;
        chgcut = big;
        dplcut = big;
        mpolecut = big;
    }
    repcut = 6.0;
    ctrncut = 6.0;
    ewaldcut = 7.0;
    dewaldcut = 7.0;
    usolvcut = 4.5;

    // set defaults for tapering, Hessian cutoff and neighbor buffers
    vdwtaper = 0.90;
    reptaper = 0.90;
    disptaper = 0.90;
    chgtaper = 0.65;
    dpltaper = 0.75;
    mpoletaper = 0.65;
    ctrntaper = 0.90;
    hesscut = 0.0;
    lbuffer = 2.0;
    pbuffer = 2.0;

    // set defaults for Ewald sum, tapering style and neighbor method
    use_ewald = false;
    use_dewald = false;
    truncate = false;
    use_lights = false;
    use_list = false;
    use_vlist = false;
    use_dlist = false;
    use_clist = false;
    use_mlist = false;
    use_ulist = false;
    dovlst = true;
    dodlst = true;
    doclst = true;
    domlst = true;
    doulst = true;

    // search the keywords for various cutoff parameters
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);

        // get values related to use of Ewald for electrostatics
        if (keyword == "EWALD") {
            use_ewald = true;
        }
        else if (keyword == "EWALD-CUTOFF") {
            iss >> ewaldcut;
        }

        // get values related to use of Ewald for dispersion
        else if (keyword == "DEWALD") {
            use_dewald = true;
        }
        else if (keyword == "DEWALD-CUTOFF") {
            iss >> dewaldcut;
        }

        // get values for the tapering style and neighbor method
        else if (keyword == "TRUNCATE") {
            truncate = true;
        }
        else if (keyword == "LIGHTS") {
            use_lights = true;
        }
        else if (keyword == "NEIGHBOR-LIST") {
            use_list = true;
            use_vlist = true;
            use_dlist = true;
            use_clist = true;
            use_mlist = true;
            use_ulist = true;
        }
        else if (keyword == "VDW-LIST") {
            use_list = true;
            use_vlist = true;
        }
        else if (keyword == "DISP-LIST") {
            use_list = true;
            use_dlist = true;
        }
        else if (keyword == "CHARGE-LIST") {
            use_list = true;
            use_clist = true;
        }
        else if (keyword == "MPOLE-LIST") {
            use_list = true;
            use_mlist = true;
            use_ulist = true;
        }

        // get values for the dipole solver preconditioner
        else if (keyword == "USOLVE-LIST") {
            use_list = true;
            use_ulist = true;
        }
        else if (keyword == "USOLVE-CUTOFF") {
            if (usolvcut != 0.) iss >> usolvcut;
        }
        else if (keyword == "USOLVE-DIAGONAL") {
            usolvcut = 0.;
        }

        // get cutoff for the magnitude of Hessian elements
        else if (keyword == "HESSIAN-CUTOFF") {
            iss >> hesscut;
        }
        // get the cutoff radii for potential energy functions
        else if (keyword == "CUTOFF") {
            iss >> value;
            vdwcut = value;
            repcut = value;
            dispcut = value;
            chgcut = value;
            dplcut = value;
            mpolecut = value;
            ewaldcut = value;
            dewaldcut = value;
            ctrncut = value;
        }
        else if (keyword == "VDW-CUTOFF") {
            iss >> vdwcut;
        }
        else if (keyword == "REPULS-CUTOFF") {
            iss >> repcut;
        }
        else if (keyword == "DISP-CUTOFF") {
            iss >> dispcut;
        }
        else if (keyword == "CHARGE-CUTOFF") {
            iss >> chgcut;
        }
        else if (keyword == "DIPOLE-CUTOFF") {
            iss >> dplcut;
        }
        else if (keyword == "MPOLE-CUTOFF") {
            iss >> mpolecut;
        }
        else if (keyword == "CHGTRN-CUTOFF") {
            iss >> ctrncut;
        }

        // get distance for initialization of energy switching
        else if (keyword == "TAPER") {
            iss >> value;
            vdwtaper = value;
            reptaper = value;
            disptaper = value;
            chgtaper = value;
            dpltaper = value;
            mpoletaper = value;
            ctrntaper = value;
        }
        else if (keyword == "VDW-TAPER") {
            iss >> vdwtaper;
        }
        else if (keyword == "REPULS-TAPER") {
            iss >> reptaper;
        }
        else if (keyword == "DISP-TAPER") {
            iss >> disptaper;
        }
        else if (keyword == "CHARGE-TAPER") {
            iss >> chgtaper;
        }
        else if (keyword == "DIPOLE-TAPER") {
            iss >> dpltaper;
        }
        else if (keyword == "MPOLE-TAPER") {
            iss >> mpoletaper;
        }
        else if (keyword == "CHGTRN-TAPER") {
            iss >> ctrntaper;
        }

        // get buffer width for use with pairwise neighbor lists
        else if (keyword == "LIST-BUFFER") {
            iss >> lbuffer;
        }
        else if (keyword == "USOLVE-BUFFER") {
            iss >> pbuffer;
        }
    }

    // check to see if preconditioner list should be disabled
    if (poltyp == "DIRECT")  use_ulist = false;
    if (usolvcut <= 0.)  use_ulist = false;
    if (use_list)  usolvcut = usolvcut - pbuffer;

    // apply any Ewald cutoff to dispersion and electrostatics
    if (use_ewald) {
        chgcut = ewaldcut;
        mpolecut = ewaldcut;
    }
    if (use_dewald) {
        dispcut = dewaldcut;
    }

    // convert any tapering percentages to absolute distances
    if (vdwtaper < 1.)  vdwtaper = vdwtaper * vdwcut;
    if (reptaper < 1.)  reptaper = reptaper * repcut;
    if (disptaper < 1.)  disptaper = disptaper * dispcut;
    if (chgtaper < 1.)  chgtaper = chgtaper * chgcut;
    if (dpltaper < 1.)  dpltaper = dpltaper * dplcut;
    if (mpoletaper < 1.)  mpoletaper = mpoletaper * mpolecut;
    if (ctrntaper < 1.)  ctrntaper = ctrntaper * ctrncut;

    // apply truncation cutoffs if they were requested
    if (truncate) {
        vdwtaper = big;
        reptaper = big;
        disptaper = big;
        chgtaper = big;
        dpltaper = big;
        mpoletaper = big;
        ctrntaper = big;
    }

    // set buffer region limits for pairwise neighbor lists
    lbuf2 = std::pow(0.5*lbuffer, 2);
    pbuf2 = std::pow(0.5*pbuffer, 2);
    vbuf2 = std::pow(vdwcut+lbuffer, 2);
    dbuf2 = std::pow(dispcut+lbuffer, 2);
    cbuf2 = std::pow(chgcut+lbuffer, 2);
    mbuf2 = std::pow(mpolecut+lbuffer, 2);
    ubuf2 = std::pow(usolvcut+pbuffer, 2);
    vbufx = std::pow(vdwcut+2.*lbuffer, 2);
    dbufx = std::pow(dispcut+2.*lbuffer, 2);
    cbufx = std::pow(chgcut+2.*lbuffer, 2);
    mbufx = std::pow(mpolecut+2.*lbuffer, 2);
    ubufx = std::pow(usolvcut+2.*pbuffer, 2);

    // specify maximum size for each of the neighbor lists
    maxvlst = 2500;
    if (vdwcut!=big and dispcut!=big) {
        limit = static_cast<int>(std::pow(std::sqrt(std::max(vbuf2,dbuf2)), 3)) + 100;
        maxvlst = std::min(limit,maxvlst);
    }
    else if (vdwcut != big) {
        limit = static_cast<int>(std::pow(std::sqrt(vbuf2), 3)) + 100;
        maxvlst = std::min(limit,maxvlst);
    }
    else if (dispcut != big) {
        limit = static_cast<int>(std::pow(std::sqrt(dbuf2), 3)) + 100;
        maxvlst = std::min(limit,maxvlst);
    }
    maxelst = 2500;
    if (chgcut!=big and mpolecut!=big) {
        limit = static_cast<int>(std::pow(std::sqrt(std::max(cbuf2,mbuf2)), 3)) + 100;
        maxelst = std::min(limit,maxelst);
    }
    else if (chgcut != big) {
        limit = static_cast<int>(std::pow(std::sqrt(cbuf2), 3)) + 100;
        maxelst = std::min(limit,maxelst);
    }
    else if (mpolecut != big) {
        limit = static_cast<int>(std::pow(std::sqrt(mbuf2), 3)) + 100;
        maxelst = std::min(limit,maxelst);
    }
    maxulst = 500;
    limit = static_cast<int>(std::pow(std::sqrt(ubuf2),3)) + 100;
    maxulst = std::min(limit,maxulst);

    // perform dynamic allocation of some global arrays
    if (use_vlist or use_dlist) {
        if (nvlst.size() != 0) nvlst.resize(0);
        if (vlst.size() != 0) vlst.resize(0);
        if (xvold.size() != 0) xvold.resize(0);
        if (yvold.size() != 0) yvold.resize(0);
        if (zvold.size() != 0) zvold.resize(0);
        nvlst.resize(n);
        vlst.resize(n, std::vector<int>(maxvlst));
        xvold.resize(n);
        yvold.resize(n);
        zvold.resize(n);
    }
    if (use_clist or use_mlist) {
        if (nelst.size() != 0) nelst.resize(0);
        if (elst.size() != 0) elst.resize(0);
        if (xeold.size() != 0) xeold.resize(0);
        if (yeold.size() != 0) yeold.resize(0);
        if (zeold.size() != 0) zeold.resize(0);
        nelst.resize(n);
        elst.resize(n,std::vector<int>(maxelst));
        xeold.resize(n);
        yeold.resize(n);
        zeold.resize(n);
        if (poltyp != "DIRECT") {
            if (tindex.size() != 0) tindex.resize(0);
            if (tdipdip.size() != 0) tdipdip.resize(0);
            tindex.resize(n*maxelst, std::vector<int>(2));
            tdipdip.resize(n*maxelst, std::vector<real>(6));
        }
    }
    if (use_ulist) {
        if (nulst.size() != 0) nulst.resize(0);
        if (ulst.size() != 0) ulst.resize(0);
        if (xuold.size() != 0) xuold.resize(0);
        if (yuold.size() != 0) yuold.resize(0);
        if (zuold.size() != 0) zuold.resize(0);
        nulst.resize(n);
        ulst.resize(n, std::vector<int>(maxulst));
        xuold.resize(n);
        yuold.resize(n);
        zuold.resize(n);
    }
}
}
