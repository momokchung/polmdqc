////////////////////////////////////////////////////
//                                                //
//  empole.cpp  --  atomic multipole calculation  //
//                                                //
////////////////////////////////////////////////////

// "empole" calculates the electrostatic energy and/or
// gradient due to atomic multipole interactions


#include "calcMode.h"
#include "empole.h"
#include "qcmdlimits.h"

void empole_a(calcMode calculationMode);
void empole_b(calcMode calculationMode);
void empole_c(calcMode calculationMode);
void empole_d(calcMode calculationMode);

void empole(calcMode calculationMode)
{
    // choose the method to sum over multipole interactions
    if (use_ewald) {
        if (use_mlist) {
            empole_d(calculationMode);
        }
        else {
            empole_c(calculationMode);
        }
    }
    else {
        if (use_mlist) {
            empole_b(calculationMode);
        }
        else {
            empole_a(calculationMode);
        }
    }
}


///////////////////////////////////////////////////////////
//                                                       //
//  empole_a.cpp  --  double loop multipole calculation  //
//                                                       //
///////////////////////////////////////////////////////////

// "empole3a" calculates the atomic multipole interaction energy
// and/or graident using a double loop


#include "action.h"
#include "analyz.h"
#include "atomid.h"
#include "atoms.h"
#include "bound.h"
#include "cell.h"
#include "chgpen.h"
#include "chgpot.h"
#include "chkpole.h"
#include "couple.h"
#include "energi.h"
#include "group.h"
#include "inform.h"
#include "inter.h"
#include "mathConst.h"
#include "molcul.h"
#include "mplpot.h"
#include "mpole.h"
#include "potent.h"
#include "rotpole.h"
#include "shunt.h"
#include "usage.h"

void empole_a(calcMode calculationMode)
{
    int i,j,k;
    int ii,kk;
    int ix,iy,iz;
    int kx,ky,kz;
    double e,f,fgrp;
    double xi,yi,zi;
    double xr,yr,zr;
    double r,r2,rr1,rr3;
    double rr5,rr7,rr9;
    double rr1i,rr3i,rr5i;
    double rr1k,rr3k,rr5k;
    double rr1ik,rr3ik,rr5ik;
    double rr7ik,rr9ik;
    double ci,dix,diy,diz;
    double qixx,qixy,qixz;
    double qiyy,qiyz,qizz;
    double ck,dkx,dky,dkz;
    double qkxx,qkxy,qkxz;
    double qkyy,qkyz,qkzz;
    double dir,dkr,dik,qik;
    double qix,qiy,qiz,qir;
    double qkx,qky,qkz,qkr;
    double diqk,dkqi,qiqk;
    double corei,corek;
    double vali,valk;
    double alphai,alphak;
    double term1,term2,term3;
    double term4,term5;
    double term1i,term2i,term3i;
    double term1k,term2k,term3k;
    double term1ik,term2ik,term3ik;
    double term4ik,term5ik;
    double dmpi[9],dmpk[9];
    double dmpik[9];
    std::vector<double> mscale;
    bool proceed;
    bool header,huge;
    bool usei,usek;
    std::string mode;
    bool do_a = false;
    bool do_g = false;
    if (calculationMode == calcMode::ANALYSIS) do_a = true;
    if (calculationMode == calcMode::GRADIENT) do_g = true;

    // zero out total atomic multipole energy and partitioning
    em = 0.;
    if (do_a) {
        nem = 0;
        for (int i = 0; i < n; i++) {
            aem[i] = 0.;
        }
    }
    if (npole == 0)  return

    // check the sign of multipole components at chiral sites
    chkpole();

    // rotate the multipole components into the global frame
    rotpole("MPOLE");
    for (int i = 0; i < n; i++) {
        printf(" %11d %20.16f\n", i+1, rpole[i][0]);
    }
// c
// c     perform dynamic allocation of some local arrays
// c
//       allocate (mscale(n))
// c
// c     initialize connected atom exclusion coefficients
// c
//       do i = 1, n
//          mscale(i) = 1.
//       end do
// c
// c     set conversion factor, cutoff and switching coefficients
// c
//       f = electric / dielec
//       mode = 'MPOLE'
//       call switch (mode)
// c
// c     print header information if debug output was requested
// c
//       header = .true.
//       if (debug .and. npole.ne.0) then
//          header = .false.
//          write (iout,10)
//    10    format (/,' Individual Atomic Multipole Interactions :',
//      &           //,' Type',14x,'Atom Names',15x,'Distance',
//      &              8x,'Energy',/)
//       end if
}



// c
// c

// c
// c     calculate the multipole interaction energy term
// c
//       do ii = 1, npole-1
//          i = ipole(ii)
//          iz = zaxis(i)
//          ix = xaxis(i)
//          iy = abs(yaxis(i))
//          xi = x(i)
//          yi = y(i)
//          zi = z(i)
//          ci = rpole(1,i)
//          dix = rpole(2,i)
//          diy = rpole(3,i)
//          diz = rpole(4,i)
//          qixx = rpole(5,i)
//          qixy = rpole(6,i)
//          qixz = rpole(7,i)
//          qiyy = rpole(9,i)
//          qiyz = rpole(10,i)
//          qizz = rpole(13,i)
//          if (use_chgpen) then
//             corei = pcore(i)
//             vali = pval(i)
//             alphai = palpha(i)
//          end if
//          usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
// c
// c     set exclusion coefficients for connected atoms
// c
//          do j = 1, n12(i)
//             mscale(i12(j,i)) = m2scale
//          end do
//          do j = 1, n13(i)
//             mscale(i13(j,i)) = m3scale
//          end do
//          do j = 1, n14(i)
//             mscale(i14(j,i)) = m4scale
//          end do
//          do j = 1, n15(i)
//             mscale(i15(j,i)) = m5scale
//          end do
// c
// c     evaluate all sites within the cutoff distance
// c
//          do kk = ii+1, npole
//             k = ipole(kk)
//             kz = zaxis(k)
//             kx = xaxis(k)
//             ky = abs(yaxis(k))
//             usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
//             proceed = .true.
//             if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
//             if (.not. use_intra)  proceed = .true.
//             if (proceed)  proceed = (usei .or. usek)
//             if (proceed) then
//                xr = x(k) - xi
//                yr = y(k) - yi
//                zr = z(k) - zi
//                if (use_bounds)  call image (xr,yr,zr)
//                r2 = xr*xr + yr* yr + zr*zr
//                if (r2 .le. off2) then
//                   r = sqrt(r2)
//                   ck = rpole(1,k)
//                   dkx = rpole(2,k)
//                   dky = rpole(3,k)
//                   dkz = rpole(4,k)
//                   qkxx = rpole(5,k)
//                   qkxy = rpole(6,k)
//                   qkxz = rpole(7,k)
//                   qkyy = rpole(9,k)
//                   qkyz = rpole(10,k)
//                   qkzz = rpole(13,k)
// c
// c     intermediates involving moments and separation distance
// c
//                   dir = dix*xr + diy*yr + diz*zr
//                   qix = qixx*xr + qixy*yr + qixz*zr
//                   qiy = qixy*xr + qiyy*yr + qiyz*zr
//                   qiz = qixz*xr + qiyz*yr + qizz*zr
//                   qir = qix*xr + qiy*yr + qiz*zr
//                   dkr = dkx*xr + dky*yr + dkz*zr
//                   qkx = qkxx*xr + qkxy*yr + qkxz*zr
//                   qky = qkxy*xr + qkyy*yr + qkyz*zr
//                   qkz = qkxz*xr + qkyz*yr + qkzz*zr
//                   qkr = qkx*xr + qky*yr + qkz*zr
//                   dik = dix*dkx + diy*dky + diz*dkz
//                   qik = qix*qkx + qiy*qky + qiz*qkz
//                   diqk = dix*qkx + diy*qky + diz*qkz
//                   dkqi = dkx*qix + dky*qiy + dkz*qiz
//                   qiqk = 2.*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
//      &                      + qixx*qkxx + qiyy*qkyy + qizz*qkzz
// c
// c     get reciprocal distance terms for this interaction
// c
//                   rr1 = f * mscale(k) / r
//                   rr3 = rr1 / r2
//                   rr5 = 3. * rr3 / r2
//                   rr7 = 5. * rr5 / r2
//                   rr9 = 7. * rr7 / r2
// c
// c     find damped multipole intermediates and energy value
// c
//                   if (use_chgpen) then
//                      corek = pcore(k)
//                      valk = pval(k)
//                      alphak = palpha(k)
//                      term1 = corei*corek
//                      term1i = corek*vali
//                      term2i = corek*dir
//                      term3i = corek*qir
//                      term1k = corei*valk
//                      term2k = -corei*dkr
//                      term3k = corei*qkr
//                      term1ik = vali*valk
//                      term2ik = valk*dir - vali*dkr + dik
//                      term3ik = vali*qkr + valk*qir - dir*dkr
//      &                            + 2.*(dkqi-diqk+qiqk)
//                      term4ik = dir*qkr - dkr*qir - 4.*qik
//                      term5ik = qir*qkr
//                      call damppole (r,9,alphai,alphak,
//      &                               dmpi,dmpk,dmpik)
//                      rr1i = dmpi(1)*rr1
//                      rr3i = dmpi(3)*rr3
//                      rr5i = dmpi(5)*rr5
//                      rr1k = dmpk(1)*rr1
//                      rr3k = dmpk(3)*rr3
//                      rr5k = dmpk(5)*rr5
//                      rr1ik = dmpik(1)*rr1
//                      rr3ik = dmpik(3)*rr3
//                      rr5ik = dmpik(5)*rr5
//                      rr7ik = dmpik(7)*rr7
//                      rr9ik = dmpik(9)*rr9
//                      e = term1*rr1 + term1i*rr1i
//      &                      + term1k*rr1k + term1ik*rr1ik
//      &                      + term2i*rr3i + term2k*rr3k
//      &                      + term2ik*rr3ik + term3i*rr5i
//      &                      + term3k*rr5k + term3ik*rr5ik
//      &                      + term4ik*rr7ik + term5ik*rr9ik
// c
// c     find standard multipole intermediates and energy value
// c
//                   else
//                      term1 = ci*ck
//                      term2 = ck*dir - ci*dkr + dik
//                      term3 = ci*qkr + ck*qir - dir*dkr
//      &                          + 2.*(dkqi-diqk+qiqk)
//                      term4 = dir*qkr - dkr*qir - 4.*qik
//                      term5 = qir*qkr
//                      e = term1*rr1 + term2*rr3 + term3*rr5
//      &                      + term4*rr7 + term5*rr9
//                   end if
// c
// c     increment the overall multipole energy components
// c
//                   if (use_group)  e = e * fgrp
//                   if (e .ne. 0.) then
//                      nem = nem + 1
//                      em = em + e
//                      aem(i) = aem(i) + 0.5*e
//                      aem(k) = aem(k) + 0.5*e
//                      if (molcule(i) .ne. molcule(k)) then
//                         einter = einter + e
//                      end if
//                   end if
// c
// c     print message if the energy of this interaction is large
// c
//                   huge = (abs(e) .gt. 100.)
//                   if ((debug.and.e.ne.0.)
//      &                  .or. (verbose.and.huge)) then
//                      if (header) then
//                         header = .false.
//                         write (iout,20)
//    20                   format (/,' Individual Atomic Multipole',
//      &                             ' Interactions :',
//      &                          //,' Type',14x,'Atom Names',
//      &                             15x,'Distance',8x,'Energy',/)
//                      end if
//                      write (iout,30)  i,name(i),k,name(k),r,e
//    30                format (' Mpole',5x,2(i7,'-',a3),9x,
//      &                          f10.4,2x,f12.4)
//                   end if
//                end if
//             end if
//          end do
// c
// c     reset exclusion coefficients for connected atoms
// c
//          do j = 1, n12(i)
//             mscale(i12(j,i)) = 1.
//          end do
//          do j = 1, n13(i)
//             mscale(i13(j,i)) = 1.
//          end do
//          do j = 1, n14(i)
//             mscale(i14(j,i)) = 1.
//          end do
//          do j = 1, n15(i)
//             mscale(i15(j,i)) = 1.
//          end do
//       end do
void empole_b(calcMode calculationMode){}
void empole_c(calcMode calculationMode){}
void empole_d(calcMode calculationMode){}
