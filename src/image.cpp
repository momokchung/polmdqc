/////////////////////////////////////////////////////////
//                                                     //
//  image.cpp  --  compute the minimum image distance  //
//                                                     //
/////////////////////////////////////////////////////////

// "image" takes the components of pairwise distance between
// two points in a periodic box and converts to the components
// of the minimum image distance
// 
// literature reference:
// 
// U. K. Deiters, "Efficient Coding of the Minimum Image Convention",
// Zeitschrift fur Physikalische Chemie, 227, 345-352 (2013)
// 
// note the "do while" clause below can be written using the "nint"
// intrinsic, and the two forms give equivalent values:
// 
// do while (abs(xr) .gt. xbox2)
//    xr = xr - sign(xbox,xr)    vs.  xr = xr - xbox*nint(xr/xbox)
// end do
// 
// which one is faster depends upon specific machine and compiler
// combinations, and other implementations are also possible


#include "boxes.h"
#include "cell.h"
#include "image.h"
#include "mathConst.h"
#include <cmath>

void image(double& xr, double& yr, double& zr)
{
    double corr;

    // for orthogonal lattice, find the desired image directly
    if (orthogonal) {
        while (std::abs(xr) > xcell2) {
            xr = xr - std::copysign(xcell,xr);
        }
        while (std::abs(yr) > ycell2) {
            yr = yr - std::copysign(ycell,yr);
        }
        while (std::abs(zr) > zcell2) {
            zr = zr - std::copysign(zcell,zr);
        }
    }

    // for monoclinic lattice, convert x and z to fractional,
    // find desired image, then translate back to Cartesian
    else if (monoclinic) {
        zr = zr / beta_sin;
        xr = xr - zr*beta_cos;
        while (std::abs(xr) > xcell2) {
            xr = xr - std::copysign(xcell,xr);
        }
        while (std::abs(yr) > ycell2) {
            yr = yr - std::copysign(ycell,yr);
        }
        while (std::abs(zr) > zcell2) {
            zr = zr - std::copysign(zcell,zr);
        }
        xr = xr + zr*beta_cos;
        zr = zr * beta_sin;
    }

    // for triclinic lattice, convert to fractional coordinates,
    // find image, then translate fractional back to Cartesian
    else if (triclinic) {
        zr = zr / gamma_term;
        yr = (yr - zr*beta_term) / gamma_sin;
        xr = xr - yr*gamma_cos - zr*beta_cos;
        while (std::abs(xr) > xcell2) {
            xr = xr - std::copysign(xcell,xr);
        }
        while (std::abs(yr) > ycell2) {
            yr = yr - std::copysign(ycell,yr);
        }
        while (std::abs(zr) > zcell2) {
            zr = zr - std::copysign(zcell,zr);
        }
        xr = xr + yr*gamma_cos + zr*beta_cos;
        yr = yr*gamma_sin + zr*beta_term;
        zr = zr * gamma_term;
    }

    // for truncated octahedron, remove the corner pieces
    else if (octahedron) {
        while (std::abs(xr) > xbox2) {
            xr = xr - std::copysign(xbox,xr);
        }
        while (std::abs(yr) > ybox2) {
            yr = yr - std::copysign(ybox,yr);
        }
        while (std::abs(zr) > zbox2) {
            zr = zr - std::copysign(zbox,zr);
        }
        if (std::abs(xr)+std::abs(yr)+std::abs(zr) > box34) {
            xr = xr - std::copysign(xbox2,xr);
            yr = yr - std::copysign(ybox2,yr);
            zr = zr - std::copysign(zbox2,zr);
        }
    }

    // for rhombic dodecahedron, align along the x- and y-axes
    else if (dodecadron) {
        while (std::abs(xr) > xbox2) {
            xr = xr - std::copysign(xbox,xr);
        }
        while (std::abs(yr) > ybox2) {
            yr = yr - std::copysign(ybox,yr);
        }
        zr = zr - root2*zbox*std::round(zr/(zbox*root2));
        corr = xbox2 * int(std::abs(xr/xbox)+std::abs(yr/ybox)+std::abs(root2*zr/zbox));
        xr = xr - std::copysign(corr,xr);
        yr = yr - std::copysign(corr,yr);
        zr = zr - std::copysign(corr,zr)*root2;
    }
}

// TODO
// c
// c     ###########################################################
// c     ##                                                       ##
// c     ##  subroutine imagen  --  fast minimum image magnitude  ##
// c     ##                                                       ##
// c     ###########################################################
// c
// c
// c     "imagen" takes the components of pairwise distance between
// c     two points and converts to the components of the minimum
// c     image distance
// c
// c     note this is a fast version for use in computing the 3D
// c     distance during neighbor list construction
// c
// c
//       subroutine imagen (xr,yr,zr)
//       use boxes
//       use math
//       implicit none
//       real*8 xr,yr,zr
//       real*8 corr
// c
// c
// c     for orthogonal lattice, find the desired image directly
// c
//       if (orthogonal) then
//          xr = xr - xbox*nint(xr/xbox)
//          yr = yr - ybox*nint(yr/ybox)
//          zr = zr - zbox*nint(zr/zbox)
// c
// c     for monoclinic lattice, convert x and z to fractional,
// c     find desired image, then translate back to Cartesian
// c
//       else if (monoclinic) then
//          zr = zr / beta_sin
//          xr = xr - zr*beta_cos
//          xr = xr - xbox*nint(xr/xbox)
//          yr = yr - ybox*nint(yr/ybox)
//          zr = zr - zbox*nint(zr/zbox)
//          xr = xr + zr*beta_cos
//          zr = zr * beta_sin
// c
// c     for triclinic lattice, convert to fractional coordinates,
// c     find image, then translate fractional back to Cartesian
// c
//       else if (triclinic) then
//          zr = zr / gamma_term
//          yr = (yr - zr*beta_term) / gamma_sin
//          xr = xr - yr*gamma_cos - zr*beta_cos
//          xr = xr - xbox*nint(xr/xbox)
//          yr = yr - ybox*nint(yr/ybox)
//          zr = zr - zbox*nint(zr/zbox)
//          xr = xr + yr*gamma_cos + zr*beta_cos
//          yr = yr*gamma_sin + zr*beta_term
//          zr = zr * gamma_term
// c
// c     for truncated octahedron, remove the corner pieces
// c
//       else if (octahedron) then
//          xr = xr - xbox*nint(xr/xbox)
//          yr = yr - ybox*nint(yr/ybox)
//          zr = zr - zbox*nint(zr/zbox)
//          if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
//             xr = xr - sign(xbox2,xr)
//             yr = yr - sign(ybox2,yr)
//             zr = zr - sign(zbox2,zr)
//          end if
// c
// c     for rhombic dodecahedron, align along the x- and y-axes
// c
//       else if (dodecadron) then
//          xr = xr - xbox*nint(xr/xbox)
//          yr = yr - ybox*nint(yr/ybox)
//          zr = zr - root2*zbox*nint(zr/(zbox*root2))
//          corr = xbox2 * int(abs(xr/xbox)+abs(yr/ybox)
//      &                        +abs(root2*zr/zbox))
//          xr = xr - sign(corr,xr)
//          yr = yr - sign(corr,yr)
//          zr = zr - sign(corr,zr)*root2
//       end if
//       return
//       end
