// Author: Moses KJ Chung
// Year:   2024

#include "atoms.h"
#include "kmulti.h"
#include "libfunc.h"
#include "mpole.h"
#include "torque.h"

namespace polmdqc
{
inline real dotp(const real a[3], const real b[3]);
inline void crossp(real ans[3], const real u[3], const real v[3]);

///////////////////////////////////////////////////////
//                                                   //
//  torque  --  convert single site torque to force  //
//                                                   //
///////////////////////////////////////////////////////

// "torque" takes the torque values on a single site defined by
// a local coordinate frame and converts to Cartesian forces on
// the original site and sites specifying the local frame, also
// gives the x,y,z-force components needed for virial computation
// 
// force distribution for the 3-fold local frame by Chao Lu,
// Ponder Lab, Washington University, July 2016
// 
// literature reference:
// 
// P. L. Popelier and A. J. Stone, "Formulae for the First and
// Second Derivatives of Anisotropic Potentials with Respect to
// Geometrical Parameters", Molecular Physics, 82, 411-425 (1994)
// 
// C. Segui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
// Representation of Electrostatics in Classical Force Fields:
// Efficient Implementation of Multipolar Interactions in
// Biomolecular Simulations", Journal of Chemical Physics, 120,
// 73-87 (2004)

template <CalcMode CalculationMode>
void torque(const std::vector<std::vector<real>>* trqPtr, std::vector<std::vector<real>>* dePtr)
{
    // integer i,j
    int ia,ib,ic,id;
    real du,dv,dw,dot;
    real usiz,vsiz,wsiz;
    real psiz,rsiz,ssiz;
    real t1siz,t2siz;
    real uvsiz,uwsiz,vwsiz;
    real ursiz,ussiz;
    real vssiz,wssiz;
    real delsiz,dphiddel;
    real uvcos,uwcos,urcos;
    real vwcos,vscos,wscos;
    real upcos,vpcos,wpcos;
    real rwcos,rucos,rvcos;
    real ut1cos,ut2cos;
    real uvsin,uwsin,ursin;
    real vwsin,vssin,wssin;
    real rwsin,rusin,rvsin;
    real ut1sin,ut2sin;
    real dphidu,dphidv,dphidw;
    real dphidr,dphids;
    real u[3],v[3],w[3];
    real p[3],r[3],s[3];
    real t1[3],t2[3];
    real uv[3],uw[3],vw[3];
    real ur[3],us[3];
    real vs[3],ws[3];
    real del[3],eps[3];
    LocalFrame axetyp;

    // choose calculation mode
    constexpr CalcFlag flags = getCalculationFlags<CalculationMode>();
    constexpr bool do_g = flags.do_gradient;
    constexpr bool do_v = flags.do_virial;

    // dereference pointers
    const std::vector<std::vector<real>>& trq = *trqPtr;
    std::vector<std::vector<real>>& de = *dePtr;

    // resolve site torques then increment forces and virial
    for (int i = 0; i < n; i++) {
        // copy trq value
        const real trqi[3] = {trq[i][0],trq[i][1],trq[i][2]};

        // initialize forces
        real frcx[3] = {0,0,0};
        real frcy[3] = {0,0,0};
        real frcz[3] = {0,0,0};

        // get the local frame type and the frame-defining atoms
        axetyp = polaxe[i];
        if (axetyp == LocalFrame::None) return;
        ia = zaxis[i] - 1;
        ib = i;
        ic = xaxis[i] - 1;
        id = std::abs(yaxis[i]) - 1;

        // construct the three rotation axes for the local frame
        u[0] = x[ia] - x[ib];
        u[1] = y[ia] - y[ib];
        u[2] = z[ia] - z[ib];
        usiz = REAL_SQRT(dotp(u,u));
        if (axetyp != LocalFrame::ZOnly) {
            v[0] = x[ic] - x[ib];
            v[1] = y[ic] - y[ib];
            v[2] = z[ic] - z[ib];
            vsiz = REAL_SQRT(dotp(v,v));
        }
        else {
            v[0] = 1;
            v[1] = 0;
            v[2] = 0;
            vsiz = 1;
            dot = u[0] / usiz;
            if (REAL_ABS(dot) > (real)0.866) {
                v[0] = 0;
                v[1] = 1;
            }
        }
        if (axetyp==LocalFrame::ZBisect or axetyp==LocalFrame::ThreeFold) {
            w[1] = x[id] - x[ib];
            w[2] = y[id] - y[ib];
            w[3] = z[id] - z[ib];
        }
        else {
            crossp(w,u,v);
        }
        wsiz = REAL_SQRT(dotp(w,w));
        for (int j = 0; j < 3; j++) {
            u[j] = u[j] / usiz;
            v[j] = v[j] / vsiz;
            w[j] = w[j] / wsiz;
        }
        // build some additional axes for the Z-Bisect local frame
        if (axetyp == LocalFrame::ZBisect) {
            r[0] = v[0] + w[0];
            r[1] = v[1] + w[1];
            r[2] = v[2] + w[2];
            rsiz = REAL_SQRT(dotp(r,r));
            crossp(s,u,r);
            ssiz = REAL_SQRT(dotp(s,s));
            for (int j = 0; j < 3; j++) {
                r[j] = r[j] / rsiz;
                s[j] = s[j] / ssiz;
            }
        }

        // negative of dot product of torque with unit vectors gives
        // result of infinitesimal rotation around these vectors
        dphidu = -dotp(trqi,u);
        dphidv = -dotp(trqi,v);
        dphidw = -dotp(trqi,w);
        if (axetyp == LocalFrame::ZBisect) {
            dphidr = -dotp(trqi,r);
            dphids = -dotp(trqi,s);
        }

        // find the perpendicular and angle for each pair of axes
        crossp(uv,v,u);
        uvsiz = REAL_SQRT(dotp(uv,uv));
        crossp(uw,w,u);
        uwsiz = REAL_SQRT(dotp(uw,uw));
        crossp(vw,w,v);
        vwsiz = REAL_SQRT(dotp(vw,vw));
        for (int j = 0; j < 3; j++) {
            uv[j] = uv[j] / uvsiz;
            uw[j] = uw[j] / uwsiz;
            vw[j] = vw[j] / vwsiz;
        }
        if (axetyp == LocalFrame::ZBisect) {
            crossp(ur,r,u);
            ursiz = REAL_SQRT(dotp(ur,ur));
            crossp(us,s,u);
            ussiz = REAL_SQRT(dotp(us,us));
            crossp(vs,s,v);
            vssiz = REAL_SQRT(dotp(vs,vs));
            crossp(ws,s,w);
            wssiz = REAL_SQRT(dotp(ws,ws));
            for (int j = 0; j < 3; j++) {
                ur[j] = ur[j] / ursiz;
                us[j] = us[j] / ussiz;
                vs[j] = vs[j] / vssiz;
                ws[j] = ws[j] / wssiz;
            }
        }

        // find sine and cosine of angles between the rotation axes
        uvcos = dotp(u,v);
        uvsin = REAL_SQRT(1 - uvcos*uvcos);
        uwcos = dotp(u,w);
        uwsin = REAL_SQRT(1 - uwcos*uwcos);
        vwcos = dotp(v,w);
        vwsin = REAL_SQRT(1 - vwcos*vwcos);
        if (axetyp == LocalFrame::ZBisect) {
            urcos = dotp(u,r);
            ursin = REAL_SQRT(1 - urcos*urcos);
            vscos = dotp(v,s);
            vssin = REAL_SQRT(1 - vscos*vscos);
            wscos = dotp(w,s);
            wssin = REAL_SQRT(1 - wscos*wscos);
        }

        // get projection of v and w onto the ru-plane for Z-Bisect
        if (axetyp == LocalFrame::ZBisect) {
            for (int j = 0; j < 3; j++) {
                t1[j] = v[j] - s[j]*vscos;
                t2[j] = w[j] - s[j]*wscos;
            }
            t1siz = REAL_SQRT(dotp(t1,t1));
            t2siz = REAL_SQRT(dotp(t2,t2));
            for (int j = 0; j < 3; j++) {
                t1[j] = t1[j] / t1siz;
                t2[j] = t2[j] / t2siz;
            }
            ut1cos = dotp(u,t1);
            ut1sin = REAL_SQRT(1 - ut1cos*ut1cos);
            ut2cos = dotp(u,t2);
            ut2sin = REAL_SQRT(1 - ut2cos*ut2cos);
        }

        // force distribution for Z-Only local coordinate frame
        if (axetyp == LocalFrame::ZOnly) {
            for (int j = 0; j < 3; j++) {
                du = uv[j]*dphidv/(usiz*uvsin) + uw[j]*dphidw/usiz;
                de[ia][j] = de[ia][j] + du;
                de[ib][j] = de[ib][j] - du;
                frcz[j] = frcz[j] + du;
            }
        }
        // force distribution for Z-then-X local coordinate frame
        else if (axetyp == LocalFrame::ZthenX) {
            for (int j = 0; j < 3; j++) {
                du = uv[j]*dphidv/(usiz*uvsin) + uw[j]*dphidw/usiz;
                dv = -uv[j]*dphidu/(vsiz*uvsin);
                de[ia][j] = de[ia][j] + du;
                de[ic][j] = de[ic][j] + dv;
                de[ib][j] = de[ib][j] - du - dv;
                frcz[j] = frcz[j] + du;
                frcx[j] = frcx[j] + dv;
            }
        }
        // force distribution for Bisector local coordinate frame
        else if (axetyp == LocalFrame::Bisector) {
            for (int j = 0; j < 3; j++) {
                du = uv[j]*dphidv/(usiz*uvsin) + (real)0.5*uw[j]*dphidw/usiz;
                dv = -uv[j]*dphidu/(vsiz*uvsin) + (real)0.5*vw[j]*dphidw/vsiz;
                de[ia][j] = de[ia][j] + du;
                de[ic][j] = de[ic][j] + dv;
                de[ib][j] = de[ib][j] - du - dv;
                frcz[j] = frcz[j] + du;
                frcx[j] = frcx[j] + dv;
            }
        }
        // force distribution for Z-Bisect local coordinate frame
        else if (axetyp == LocalFrame::ZBisect) {
            for (int j = 0; j < 3; j++) {
                du = ur[j]*dphidr/(usiz*ursin) + us[j]*dphids/usiz;
                dv = (vssin*s[j]-vscos*t1[j])*dphidu / (vsiz*(ut1sin+ut2sin));
                dw = (wssin*s[j]-wscos*t2[j])*dphidu / (wsiz*(ut1sin+ut2sin));
                de[ia][j] = de[ia][j] + du;
                de[ic][j] = de[ic][j] + dv;
                de[id][j] = de[id][j] + dw;
                de[ib][j] = de[ib][j] - du - dv - dw;
                frcz[j] = frcz[j] + du;
                frcx[j] = frcx[j] + dv;
                frcy[j] = frcy[j] + dw;
            }
        }
        // force distribution for 3-Fold local coordinate frame
        else if (axetyp == LocalFrame::ThreeFold) {
            p[0] = u[0] + v[0] + w[0];
            p[1] = u[1] + v[1] + w[1];
            p[2] = u[2] + v[2] + w[2];
            psiz = REAL_SQRT(dotp(p,p));
            for (int j = 0; j < 3; j++) { 
                p[j] = p[j] / psiz;
            }
            wpcos = dotp(w,p);
            upcos = dotp(u,p);
            vpcos = dotp(v,p);
            r[0] = u[0] + v[0];
            r[1] = u[1] + v[1];
            r[2] = u[2] + v[2];
            rsiz = REAL_SQRT(dotp(r,r));
            for (int j = 0; j < 3; j++) {
                r[j] = r[j] / rsiz;
            }
            rwcos = dotp(r,w);
            rwsin = REAL_SQRT(1 - rwcos*rwcos);
            dphidr = -dotp(trqi,r);
            crossp(del,r,w);
            delsiz = REAL_SQRT(dotp(del,del));
            for (int j = 0; j < 3; j++) {
                del[j] = del[j] / delsiz;
            }
            dphiddel = -dotp(trqi,del);
            crossp(eps,del,w);
            for (int j = 0; j < 3; j++) {
                dw = del[j]*dphidr/(wsiz*rwsin) + eps[j]*dphiddel*wpcos/(wsiz*psiz) ;
                de[id][j] = de[id][j] + dw;
                de[ib][j] = de[ib][j] - dw;
                frcy[j] = frcy[j] + dw;
            }
            r[0] = v[0] + w[0];
            r[1] = v[1] + w[1];
            r[2] = v[2] + w[2];
            rsiz = REAL_SQRT(dotp(r,r));
            for (int j = 0; j < 3; j++) {
                r[j] = r[j] / rsiz;
            }
            rucos = dotp(r,u);
            rusin = REAL_SQRT(1 - rucos*rucos);
            dphidr = -dotp(trqi,r);
            crossp(del,r,u);
            delsiz = REAL_SQRT(dotp(del,del));
            for (int j = 0; j < 3; j++) {
                del[j] = del[j] / delsiz;
            }
            dphiddel = -dotp(trqi,del);
            crossp(eps,del,u);
            for (int j = 0; j < 3; j++) {
                du = del[j]*dphidr/(usiz*rusin) + eps[j]*dphiddel*upcos/(usiz*psiz);
                de[ia][j] = de[ia][j] + du;
                de[ib][j] = de[ib][j] - du;
                frcz[j] = frcz[j] + du;
            }
            r[0] = u[0] + w[0];
            r[1] = u[1] + w[1];
            r[2] = u[2] + w[2];
            rsiz = REAL_SQRT(dotp(r,r));
            for (int j = 0; j < 3; j++) {
                r[j] = r[j] / rsiz;
            }
            rvcos = dotp(r,v);
            rvsin = REAL_SQRT(1 - rvcos*rvcos);
            dphidr = -dotp(trqi,r);
            crossp(del,r,v);
            delsiz = REAL_SQRT(dotp(del,del));
            for (int j = 0; j < 3; j++) {
                del[j] = del[j] / delsiz;
            }
            dphiddel = -dotp(trqi,del);
            crossp(eps,del,v);
            for (int j = 0; j < 3; j++) {
                dv = del[j]*dphidr/(vsiz*rvsin) + eps[j]*dphiddel*vpcos/(vsiz*psiz);
                de[ic][j] = de[ic][j] + dv;
                de[ib][j] = de[ib][j] - dv;
                frcx[j] = frcx[j] + dv;
            }
        }
    }
}


    // // resolve site torques then increment forces and virial
    //     int iz = ia;
    //     int ix = ic;
    //     int iy = id;
    //     if (iz == -1) iz = i;
    //     if (ix == -1) ix = i;
    //     if (iy == -1) iy = i;
    //     xiz = x[iz] - x[i];
    //     yiz = y[iz] - y[i];
    //     ziz = z[iz] - z[i];
    //     xix = x[ix] - x[i];
    //     yix = y[ix] - y[i];
    //     zix = z[ix] - z[i];
    //     xiy = x[iy] - x[i];
    //     yiy = y[iy] - y[i];
    //     ziy = z[iy] - z[i];
    //     vxx = xix*frcx[0] + xiy*frcy[0] + xiz*frcz[0]
    //     vxy = (real)0.5 * (yix*frcx[0] + yiy*frcy[0] + yiz*frcz[0]
    //  &                    + xix*frcx[1] + xiy*frcy[1] + xiz*frcz[1])
    //     vxz = (real)0.5 * (zix*frcx[0] + ziy*frcy[0] + ziz*frcz[0]
    //  &                    + xix*frcx[2] + xiy*frcy[2] + xiz*frcz[2]) 
    //     vyy = yix*frcx[1] + yiy*frcy[1] + yiz*frcz[1]
    //     vyz = (real)0.5 * (zix*frcx[1] + ziy*frcy[1] + ziz*frcz[1]
    //  &                    + yix*frcx[2] + yiy*frcy[2] + yiz*frcz[2])
    //     vzz = zix*frcx[2] + ziy*frcy[2] + ziz*frcz[2]
    //     vir(1,1) = vir(1,1) + vxx
    //     vir(2,1) = vir(2,1) + vxy
    //     vir(3,1) = vir(3,1) + vxz
    //     vir(1,2) = vir(1,2) + vxy
    //     vir(2,2) = vir(2,2) + vyy
    //     vir(3,2) = vir(3,2) + vyz
    //     vir(1,3) = vir(1,3) + vxz
    //     vir(2,3) = vir(2,3) + vyz
    //     vir(3,3) = vir(3,3) + vzz
    // end do

// define some functions to be used in torque
inline real dotp(const real a[3], const real b[3])
{
   return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline void crossp(real ans[3], const real u[3], const real v[3])
{
   ans[0] = u[1] * v[2] - u[2] * v[1];
   ans[1] = u[2] * v[0] - u[0] * v[2];
   ans[2] = u[0] * v[1] - u[1] * v[0];
}

// explicit instatiation
template void torque<CalcMode::Gradient>(const std::vector<std::vector<real>>* trqPtr, std::vector<std::vector<real>>* dePtr);
template void torque<CalcMode::Virial>(const std::vector<std::vector<real>>* trqPtr, std::vector<std::vector<real>>* dePtr);
}
