// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "bath.h"
#include "couple.h"
#include "getnumb.h"
#include "gettext.h"
#include "getword.h"
#include "gkstuf.h"
#include "inform.h"
#include "keys.h"
#include "ksolut.h"
#include "ksolv.h"
#include "kvdws.h"
#include "mathConst.h"
#include "nonpol.h"
#include "polar.h"
#include "polopt.h"
#include "polpot.h"
#include "potent.h"
#include "ptable.h"
#include "sizes.h"
#include "solpot.h"
#include "solute.h"
#include "upcase.h"
#include "vdw.h"
#include "vdwpot.h"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////
//                                             //
//  ksolv  --  solvation parameter assignment  //
//                                             //
/////////////////////////////////////////////////

// "ksolv" assigns implicit solvation energy parameters for
// the surface area, generalized Born, generalized Kirkwood,
// Poisson-Boltzmann, cavity-dispersion and HPMF models

void ksolv()
{
    int i,k,next;
    double pbrd,csrd,gkrd;
    bool header;
    std::string keyword;
    std::string value;
    std::string record;
    std::string string;
    std::istringstream iss;

    // defaults for implicit solvation term and parameters
    use_solv = false;
    use_born = false;
    solvtyp = "";
    borntyp = "";
    doffset = 0.09;
    onipr = 0.;

    // search keywords for implicit solvation commands
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "SOLVATE") {
            use_solv = true;
            use_born = false;
            solvtyp = "ASP";
            getword(record,value,next);
            upcase(value);
            if (value == "ASP") {
               solvtyp = "ASP";
            }
            else if (value == "SASA") {
               solvtyp = "SASA";
            }
            else if (value == "ONION") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "ONION";
            }
            else if (value == "STILL") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "STILL";
            }
            else if (value == "HCT") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "HCT";
            }
            else if (value == "OBC") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "OBC";
            }
            else if (value == "ACE") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "ACE";
            }
            else if (value == "GB-HPMF") {
               use_born = true;
               solvtyp = "GB-HPMF";
               borntyp = "STILL";
            }
            else if (value == "GB") {
               use_born = true;
               solvtyp = "GB";
               borntyp = "STILL";
            }
            else if (value == "GK-HPMF") {
               use_born = true;
               solvtyp = "GK-HPMF";
               borntyp = "GRYCUK";
            }
            else if (value == "GK") {
               use_born = true;
               solvtyp = "GK";
               borntyp = "GRYCUK";
            }
            else if (value == "PB-HPMF") {
               solvtyp = "PB-HPMF";
            }
            else if (value == "PB") {
               solvtyp = "PB";
            }
        }
        else if (keyword == "BORN-RADIUS") {
            getword(record,value,next);
            upcase(value);
            if (value == "ONION") {
               borntyp = "ONION";
            }
            else if (value == "STILL") {
               borntyp = "STILL";
            }
            else if (value == "HCT") {
               borntyp = "HCT";
            }
            else if (value == "OBC") {
               borntyp = "OBC";
            }
            else if (value == "ACE") {
               borntyp = "ACE";
            }
            else if (value == "GRYCUK") {
               borntyp = "GRYCUK";
            }
            else if (value == "GONION") {
               borntyp = "GONION";
            }
            else if (value == "PERFECT") {
               borntyp = "PERFECT";
            }
        }
        else if (keyword == "ONION-PROBE") {
            if (!(iss >> onipr)) break;
        }
        else if (keyword == "DIELECTRIC-OFFSET") {
            if (!(iss >> doffset)) break;
            if (doffset < 0.)  doffset = -doffset;
        }
    }

    // process keywords containing solvation parameters
    header = true;
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "SOLUTE") {
            getnumb(record,k,next);
            if (k>=1 and k<=maxtyp) {
                int km = k-1;
                pbrd = 0.;
                csrd = 0.;
                gkrd = 0.;
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> pbrd >> csrd >> gkrd;
                if (header and !silent) {
                    header = false;
                    printf("\n Additional Solvation Parameters :\n\n");
                    printf("     Atom Type                PB Size     CS Size     GK Size\n\n");
               }
               pbr[km] = pbrd;
               gkr[km] = csrd;
               pbr[km] = gkrd;
               if (!silent) {
                    printf("      %6d       %12.4f%12.4f%12.4f\n", k,pbrd,csrd,gkrd);
               }
            }
            else if (k > maxtyp) {
                printf("\n KSOLV  --  Only Atom Types through%5d are Allowed\n", maxtyp);
                informAbort = true;
            }
        }
    }

    // perform dynamic allocation of some global arrays
    if (rsolv.size() != 0) rsolv.resize(0);
    rsolv.resize(n);

    // invoke the setup needed for perfect Born radius model
    if (borntyp == "PERFECT") kpb();

    // invoke the setup needed for specific solvation models
    if (solvtyp=="ASP" or solvtyp=="SASA") {
        ksa();
    }
    else if (solvtyp == "GB-HPMF") {
        kgb();
        khpmf();
    }
    else if (solvtyp == "GB") {
        kgb();
    }
    else if (solvtyp == "GK-HPMF") {
        kgk();
        khpmf();
    }
    else if (solvtyp == "GK") {
        kgk();
        // knp();
    }
    else if (solvtyp == "PB-HPMF") {
        kpb();
        khpmf();
    }
    else if (solvtyp == "PB") {
        kpb();
        knp();
    }
}

//////////////////////////////////////////////////////
//                                                  //
//  ksa  --  set surface area solvation parameters  //
//                                                  //
//////////////////////////////////////////////////////

// "ksa" initializes parameters needed for surface area-based
// implicit solvation models including ASP and SASA
// 
// literature references:
// 
// L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
// Applied to Molecular Dynamics of Proteins in Solution",
// Protein Science, 1, 227-235 (1992)  (Eisenberg-McLachlan ASP)
// 
// T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, "Accessible
// Surface Areas as a Measure of the Thermodynamic Parameters of
// Hydration of Peptides", PNAS, 84, 3086-3090 (1987)  (SASA)

void ksa()
{

}
//       use sizes
//       use atomid
//       use atoms
//       use couple
//       use solpot
//       use solute
//       implicit none
//       integer i,j,k
//       integer atmnum
// c
// c
// c     perform dynamic allocation of some global arrays
// c
//       if (allocated(asolv))  deallocate (asolv)
//       allocate (asolv(n))
// c
// c     assign the Eisenberg-McLachlan ASP solvation parameters;
// c     parameters only available for protein-peptide groups
// c
//       if (solvtyp == "ASP") {
//          do i = 1, n
//             atmnum = atomic[i]
//             if (atmnum == 6) {
//                rsolv[i] = 1.9
//                asolv[i] = 0.004
//             else if (atmnum == 7) {
//                rsolv[i] = 1.7
//                asolv[i] = -0.113
//                if (n12[i] == 4) {
//                   asolv[i] = -0.169
//                }
//             else if (atmnum == 8) {
//                rsolv[i] = 1.4
//                asolv[i] = -0.113
//                if (n12[i]==1 and atomic(i12(1,i))==6) {
//                   do j = 1, n13[i]
//                      k = i13(j,i)
//                      if (n12[k]==1 and atomic[k]==8) {
//                         asolv[i] = -0.166
//                      }
//                   }
//                }
//                do j = 1, n12[i]
//                   k = i12(j,i)
//                   if (atomic[k] == 15)  asolv[i] = -0.14
//                }
//             else if (atmnum == 15) {
//                rsolv[i] = 1.9
//                asolv[i] = -0.14
//             else if (atmnum == 16) {
//                rsolv[i] = 1.8
//                asolv[i] = -0.017
//             else
//                rsolv[i] = 0.
//                asolv[i] = 0.
//             }
//          }
//       }
// c
// c     assign the Ooi-Scheraga SASA solvation parameters;
// c     parameters only available for protein-peptide groups
// c
//       if (solvtyp == "SASA") {
//          do i = 1, n
//             atmnum = atomic[i]
//             if (atmnum == 6) {
//                rsolv[i] = 2.
//                asolv[i] = 0.008
//                if (n12[i] == 3) {
//                   rsolv[i] = 1.75
//                   asolv[i] = -0.008
//                   do j = 1, n12[i]
//                      k = i12(j,i)
//                      if (atomic[k] == 8) {
//                         rsolv[i] = 1.55
//                         asolv[i] = 0.427
//                      }
//                   }
//                }
//             else if (atmnum == 7) {
//                rsolv[i] = 1.55
//                asolv[i] = -0.132
//                if (n12[i] == 4)  asolv[i] = -1.212
//             else if (atmnum == 8) {
//                rsolv[i] = 1.4
//                if (n12[i] == 1) {
//                   asolv[i] = -0.038
//                   if (atomic(i12(1,i)) == 6) {
//                      do j = 1, n13[i]
//                         k = i13(j,i)
//                         if (n12[k]==1 and atomic[k]==8) {
//                            asolv[i] = -0.77
//                         }
//                      }
//                   }
//                else if (n12[i] == 2) {
//                   asolv[i] = -0.172
//                }
//                do j = 1, n12[i]
//                   k = i12(j,i)
//                   if (atomic[k] == 15)  asolv[i] = -0.717
//                }
//             else if (atmnum == 15) {
//                rsolv[i] = 2.1
//                asolv[i] = 0.
//             else if (atmnum == 16) {
//                rsolv[i] = 2.
//                asolv[i] = -0.021
//             else if (atmnum == 17) {
//                rsolv[i] = 2.
//                asolv[i] = 0.012
//             else
//                rsolv[i] = 0.
//                asolv[i] = 0.
//             }
//          }
//       }

///////////////////////////////////////////////////
//                                               //
//  kgb  --  assign generalized Born parameters  //
//                                               //
///////////////////////////////////////////////////

// "kgb" initializes parameters needed for the generalized
// Born implicit solvation models
// 
// literature references:
// 
// M. Schaefer, C. Bartels, F. Leclerc and M. Karplus, "Effective
// Atom Volumes for Implicit Solvent Models: Comparison between
// Voronoi Volumes and Minimum Fluctuations Volumes", Journal of
// Computational Chemistry, 22, 1857-1879 (2001)  (ACE)


void kgb()
{
    
}
//       use sizes
//       use angbnd
//       use atmlst
//       use atomid
//       use atoms
//       use bndstr
//       use chgpot
//       use couple
//       use math
//       use potent
//       use solpot
//       use solute
//       implicit none
//       integer i,j,k,m
//       integer mm,nh,kc
//       integer ia,ib,ic,id
//       integer atmnum,atmmas
//       real*8 ri,ri2,rk,rk2
//       real*8 c1,c2,c3,pi2
//       real*8 r,r2,r4,rab,rbc
//       real*8 cosine,factor
//       real*8 h,ratio,term
//       real*8 width,qterm,temp
//       real*8 alpha,alpha2,alpha4
//       real*8 vk,prod2,prod4
//       real*8 fik,tik2,qik,uik
//       real*8 s2ik,s3ik,omgik
//       logical amide
// c
// c
// c     perform dynamic allocation of some global arrays
// c
//       if (.not. allocated(wace))  allocate (wace(maxclass,maxclass))
//       if (.not. allocated(s2ace))  allocate (s2ace(maxclass,maxclass))
//       if (.not. allocated(uace))  allocate (uace(maxclass,maxclass))
//       if (allocated(asolv))  deallocate (asolv)
//       if (allocated(rborn))  deallocate (rborn)
//       if (allocated(drb))  deallocate (drb)
//       if (allocated(drobc))  deallocate (drobc)
//       if (allocated(gpol))  deallocate (gpol)
//       if (allocated(shct))  deallocate (shct)
//       if (allocated(aobc))  deallocate (aobc)
//       if (allocated(bobc))  deallocate (bobc)
//       if (allocated(gobc))  deallocate (gobc)
//       if (allocated(vsolv))  deallocate (vsolv)
//       allocate (asolv(n))
//       allocate (rborn(n))
//       allocate (drb(n))
//       allocate (drobc(n))
//       allocate (gpol(n))
//       allocate (shct(n))
//       allocate (aobc(n))
//       allocate (bobc(n))
//       allocate (gobc(n))
//       allocate (vsolv(n))
// c
// c     set offset and scaling values for analytical Still method
// c
//       if (borntyp == "STILL") {
//          p1 = 0.073
//          p2 = 0.921
//          p3 = 6.211
//          p4 = 15.236
//          p5 = 1.254
//          if (!use_bond)  call kbond
//          if (!use_angle)  call kangle
//       }
// c
// c     set overlap scale factors for HCT and OBC methods
// c
//       if (borntyp=="HCT" or borntyp=="OBC") {
//          do i = 1, n
//             shct[i] = 0.8
//             atmnum = atomic[i]
//             if (atmnum == 1)  shct[i] = 0.85
//             if (atmnum == 6)  shct[i] = 0.72
//             if (atmnum == 7)  shct[i] = 0.79
//             if (atmnum == 8)  shct[i] = 0.85
//             if (atmnum == 9)  shct[i] = 0.88
//             if (atmnum == 15)  shct[i] = 0.86
//             if (atmnum == 16)  shct[i] = 0.96
//             if (atmnum == 26)  shct[i] = 0.88
//          }
//       }
// c
// c     set rescaling coefficients for the OBC method
// c
//       if (borntyp == "OBC") {
//          do i = 1, n
//             aobc[i] = 1.0
//             bobc[i] = 0.8
//             gobc[i] = 4.85
//          }
//       }
// c
// c     set the Gaussian width factor for the ACE method
// c
//       if (borntyp == "ACE") {
//          width = 1.2
//       }
// c
// c     assign surface area factors for nonpolar solvation
// c
//       if (borntyp == "ONION") {
//          do i = 1, n
//             asolv[i] = 0.0072
//          }
//       else if (borntyp == "STILL") {
//          do i = 1, n
//             asolv[i] = 0.0049
//          }
//       else if (borntyp == "HCT") {
//          do i = 1, n
//             asolv[i] = 0.0054
//          }
//       else if (borntyp == "OBC") {
//          do i = 1, n
//             asolv[i] = 0.0054
//          }
//       else if (borntyp == "ACE") {
//          do i = 1, n
//             asolv[i] = 0.003
//          }
//       }
// c
// c     assign standard radii for GB/SA methods other than ACE;
// c     taken from Macromodel and OPLS-AA, except for hydrogens
// c
//       if (borntyp != "ACE") {
//          do i = 1, n
//             atmnum = atomic[i]
//             if (atmnum == 1) {
//                rsolv[i] = 1.25
//                k = i12(1,i)
//                if (atomic[k] == 7)  rsolv[i] = 1.15
//                if (atomic[k] == 8)  rsolv[i] = 1.05
//             else if (atmnum == 3) {
//                rsolv[i] = 1.432
//             else if (atmnum == 6) {
//                rsolv[i] = 1.9
//                if (n12[i] == 3)  rsolv[i] = 1.875
//                if (n12[i] == 2)  rsolv[i] = 1.825
//             else if (atmnum == 7) {
//                rsolv[i] = 1.7063
//                if (n12[i] == 4)  rsolv[i] = 1.625
//                if (n12[i] == 1)  rsolv[i] = 1.6
//             else if (atmnum == 8) {
//                rsolv[i] = 1.535
//                if (n12[i] == 1)  rsolv[i] = 1.48
//             else if (atmnum == 9) {
//                rsolv[i] = 1.47
//             else if (atmnum == 10) {
//                rsolv[i] = 1.39
//             else if (atmnum == 11) {
//                rsolv[i] = 1.992
//             else if (atmnum == 12) {
//                rsolv[i] = 1.7
//             else if (atmnum == 14) {
//                rsolv[i] = 1.8
//             else if (atmnum == 15) {
//                rsolv[i] = 1.87
//             else if (atmnum == 16) {
//                rsolv[i] = 1.775
//             else if (atmnum == 17) {
//                rsolv[i] = 1.735
//             else if (atmnum == 18) {
//                rsolv[i] = 1.7
//             else if (atmnum == 19) {
//                rsolv[i] = 2.123
//             else if (atmnum == 20) {
//                rsolv[i] = 1.817
//             else if (atmnum == 35) {
//                rsolv[i] = 1.9
//             else if (atmnum == 36) {
//                rsolv[i] = 1.812
//             else if (atmnum == 37) {
//                rsolv[i] = 2.26
//             else if (atmnum == 53) {
//                rsolv[i] = 2.1
//             else if (atmnum == 54) {
//                rsolv[i] = 1.967
//             else if (atmnum == 55) {
//                rsolv[i] = 2.507
//             else if (atmnum == 56) {
//                rsolv[i] = 2.188
//             else
//                rsolv[i] = 2.
//             }
//          }
//       }
// c
// c     compute the atomic volumes for the analytical Still method
// c
//       if (borntyp == "STILL") {
//          do i = 1, n
//             vsolv[i] = (4.*pi/3.) * std::pow(rsolv[i],3)
//             ri = rsolv[i]
//             ri2 = ri * ri
//             do j = 1, n12[i]
//                k = i12(j,i)
//                rk = rsolv[k]
//                r = 1.01 * bl(bndlist(j,i))
//                ratio = (rk*rk-ri2-r*r) / (2.*ri*r)
//                h = ri * (1.+ratio)
//                term = (pi/3.) * h * h * (3.*ri-h)
//                vsolv[i] = vsolv[i] - term
//             }
//          }
// c
// c     get self-, 1-2 and 1-3 polarization for analytical Still method
// c
//          do i = 1, n
//             gpol[i] = -0.5 * electric / (rsolv[i]-doffset+p1)
//          }
//          do i = 1, nbond
//             ia = ibnd(1,i)
//             ib = ibnd(2,i)
//             r = bl[i]
//             r4 = std::pow(r,4)
//             gpol(ia) = gpol(ia) + p2*vsolv(ib)/r4
//             gpol(ib) = gpol(ib) + p2*vsolv(ia)/r4
//          }
//          do i = 1, nangle
//             ia = iang(1,i)
//             ib = iang(2,i)
//             ic = iang(3,i)
//             factor = 1.
//             do j = 1, n12(ia)
//                id = i12(j,ia)
//                if (id == ic) {
//                   factor = 0.
//                else if (id != ib) {
//                   do k = 1, n12(ic)
//                      if (i12(k,ic) == id) {
//                         factor = 0.5
//                      }
//                   }
//                }
//             }
//             do j = 1, n12(ib)
//                if (i12(j,ib) == ia) {
//                   rab = bl(bndlist(j,ib))
//                else if (i12(j,ib) == ic) {
//                   rbc = bl(bndlist(j,ib))
//                }
//             }
//             cosine = cos(anat[i]/radian)
//             r2 = std::pow(rab,2) + std::pow(rbc,2) - 2.*rab*rbc*cosine
//             r4 = r2 * r2
//             gpol(ia) = gpol(ia) + factor*p3*vsolv(ic)/r4
//             gpol(ic) = gpol(ic) + factor*p3*vsolv(ia)/r4
//          }
//       }
// c
// c     assign the atomic radii and volumes for the ACE method;
// c     volumes taken from average Voronoi values with hydrogens
// c
//       if (borntyp == "ACE") {
//          do i = 1, n
//             atmnum = atomic[i]
//             atmmas = nint(mass[i])
//             if (atmnum == 1) {
//                rsolv[i] = 1.468
//                vsolv[i] = 11.
//                k = i12(1,i)
//                if (atomic[k]==6 and n12[k]==4) {
//                   vsolv[i] = 11.895
//                else if (atomic[k]==6 and n12[k]==3) {
//                   vsolv[i] = 13.242
//                else if (atomic[k]==7 and n12[k]==4) {
//                   rsolv[i] = 0.6
//                   vsolv[i] = 9.138
//                else if (atomic[k]==7 or atomic[k]==8) {
//                   rsolv[i] = 0.6
//                   vsolv[i] = 9.901
//                else if (atomic[k]!=16) {
//                   rsolv[i] = 1.468
//                   vsolv[i] = 13.071
//                }
//             else if (atmnum == 6) {
//                rsolv[i] = 2.49
//                vsolv[i] = 7.
//                nh = 0
//                do j = 1, n12[i]
//                   k = i12(j,i)
//                   if (atomic[k] == 1)  nh = nh + 1
//                }
//                if (n12[i] == 4) {
//                   if (nh == 3) {
//                      vsolv[i] = 3.042
//                   else if (nh == 2) {
//                      vsolv[i] = 3.743
//                   else if (nh == 1) {
//                      vsolv[i] = 4.38
//                   }
//                else if (n12[i] == 3) {
//                   if (nh == 1) {
//                      rsolv[i] = 2.1
//                      vsolv[i] = 7.482
//                   else if (nh == 0) {
//                      rsolv[i] = 2.1
//                      vsolv[i] = 8.288
//                   }
//                   do j = 1, n12[i]
//                      k = i12(1,j)
//                      if (atomic[k]==8 and n12[k]==1) {
//                         rsolv[i] = 2.1
//                         vsolv[i] = 7.139
//                      }
//                   }
//                }
//                if (atmmas == 15) {
//                   rsolv[i] = 2.165
//                   vsolv[i] = 33.175
//                else if (atmmas == 14) {
//                   rsolv[i] = 2.235
//                   vsolv[i] = 20.862
//                else if (atmmas==13 and n12[i]==2) {
//                   rsolv[i] = 2.1
//                   vsolv[i] = 20.329
//                else if (atmmas == 13) {
//                   rsolv[i] = 2.365
//                   vsolv[i] = 11.784
//                }
//             else if (atmnum == 7) {
//                rsolv[i] = 1.6
//                vsolv[i] = 6.
//                nh = 0
//                do j = 1, n12[i]
//                   k = i12(j,i)
//                   if (atomic[k] == 1)  nh = nh + 1
//                }
//                if (n12[i] == 4) {
//                   if (nh == 3) {
//                      vsolv[i] = 2.549
//                   else if (nh == 2) {
//                      vsolv[i] = 3.304
//                   }
//                else if (n12[i] == 3) {
//                   amide = false
//                   do j = 1, n12[i]
//                      m = i12(j,i)
//                      if (atomic(m) == 6) {
//                         do k = 1, n12(m)
//                            mm = i12(k,m)
//                            if (atomic(mm)==8 and n12(mm)==1) {
//                               amide = true
//                            }
//                         }
//                      }
//                   }
//                   if (amide) {
//                      if (nh == 0) {
//                         vsolv[i] = 7.189
//                      else if (nh == 1) {
//                         vsolv[i] = 6.03
//                      else if (nh == 2) {
//                         vsolv[i] = 5.693
//                      }
//                   else
//                      if (nh == 2) {
//                         vsolv[i] = 5.677
//                      else if (nh == 2) {
//                         vsolv[i] = 6.498
//                      }
//                   }
//                }
//             else if (atmnum == 8) {
//                rsolv[i] = 1.6
//                vsolv[i] = 12.
//                if (n12[i] == 1) {
//                   vsolv[i] = 13.532
//                   k = i12(1,i)
//                   if (atomic[k] == 15) {
//                      vsolv[i] = 17.202
//                   else
//                      do j = 1, n13[i]
//                         k = i13(j,i)
//                         if (atomic(j)==8 and n12(j)==1) {
//                            vsolv[i] = 15.40
//                         }
//                      }
//                   }
//                else if (n12[i] == 2) {
//                   vsolv[i] = 10.642
//                   do j = 1, n12[i]
//                      k = i12(j,i)
//                      if (atomic[k] == 15)  vsolv[i] = 11.416
//                   }
//                }
//             else if (atmnum == 12) {
//                rsolv[i] = 1.
//                vsolv[i] = 15.235
//             else if (atmnum == 15) {
//                rsolv[i] = 1.89
//                vsolv[i] = 6.131
//             else if (atmnum == 16) {
//                rsolv[i] = 1.89
//                vsolv[i] = 17.232
//                do j = 1, n12[i]
//                   k = i12(j,i)
//                   if (atomic[k] == 16)  vsolv[i] = 18.465
//                }
//             else if (atmnum == 26) {
//                rsolv[i] = 0.65
//                vsolv[i] = 9.951
//             else
//                rsolv[i] = 0.
//                vsolv[i] = 0.
//             }
//          }
// c
// c     calculate the pairwise parameters for the ACE method
// c
//          c1 = 4. / (3.*pi)
//          c2 = 77. * pi * root2 / 512.
//          c3 = 2. * pi * rootpi
//          pi2 = 1. / (pi*pi)
//          do i = 1, n
//             ic = class[i]
//             ri = rsolv[i]
//             ri2 = ri * ri
//             do k = 1, n
//                kc = class[k]
//                rk = rsolv[k]
//                vk = vsolv(kc)
//                rk2 = rk * rk
//                alpha = max(width,ri/rk)
//                alpha2 = alpha * alpha
//                alpha4 = alpha2 * alpha2
//                prod2 = alpha2 * rk2
//                prod4 = prod2 * prod2
//                ratio = alpha2 * rk2 / ri2
//                tik2 = 0.5 * pi * ratio
//                temp = 1. / (1.+2.*tik2)
//                fik = 2./(1.+tik2) - temp
//                qik = tik2 * std::sqrt(temp)
//                qterm = qik - atan(qik)
//                if (k != i) {
//                   omgik = vk * qterm * pi2 / prod4
//                else
//                   omgik = c1 * qterm / (alpha4 * ri)
//                }
//                s2ik = 3. * qterm * prod2
//      &                   / ((3.+fik)*qik-4.*atan(qik))
//                s3ik = s2ik * std::sqrt(s2ik)
//                uik = c2 * ri / (1.-(c3*s3ik*ri*omgik/vk))
//                wace(ic,kc) = omgik
//                s2ace(ic,kc) = s2ik
//                uace(ic,kc) = uik
//             }
//          }
//       }

////////////////////////////////////////////////////
//                                                //
//  kgk  --  set generalized Kirkwood parameters  //
//                                                //
////////////////////////////////////////////////////

// "kgk" initializes parameters needed for the generalized
// Kirkwood implicit solvation model

void kgk()
{
    int i,it,next;
    int atmnum;
    double dhct;
    bool descreenVDW;
    bool descreenHydrogen;
    std::string radtyp;
    std::string keyword;
    std::string value;
    std::string record;
    std::string string;
    std::istringstream iss;

    // perform dynamic allocation of some global arrays
    if (rsolv.size() != 0) rsolv.resize(0);
    if (rdescr.size() != 0) rdescr.resize(0);
    if (rborn.size() != 0) rborn.resize(0);
    if (drb.size() != 0) drb.resize(0);
    if (drbp.size() != 0) drbp.resize(0);
    if (drobc.size() != 0) drobc.resize(0);
    if (shct.size() != 0) shct.resize(0);
    if (udirs.size() != 0) udirs.resize(0);
    if (udirps.size() != 0) udirps.resize(0);
    if (uinds.size() != 0) uinds.resize(0);
    if (uinps.size() != 0) uinps.resize(0);
    if (uopts.size() != 0) uopts.resize(0);
    if (uoptps.size() != 0) uoptps.resize(0);
    rsolv.resize(n);
    rdescr.resize(n);
    rborn.resize(n);
    drb.resize(n);
    drbp.resize(n);
    drobc.resize(n);
    shct.resize(n);
    udirs.resize(n, std::vector<double>(3));
    udirps.resize(n, std::vector<double>(3));
    uinds.resize(n, std::vector<double>(3));
    uinps.resize(n, std::vector<double>(3));
    if (poltyp == "OPT") {
        uopts.resize(n, std::vector<std::vector<double>>(3, std::vector<double>(optorder+1)));
        uoptps.resize(n, std::vector<std::vector<double>>(3, std::vector<double>(optorder+1)));
    }

    // set default value for exponent in the GB/GK function
    gkc = 2.455;
    dhct = 0.72;
    radtyp = "SOLUTE";
    descreenVDW = true;
    descreenHydrogen = false;

    // get any altered generalized Kirkwood values from keyfile
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext (record,keyword,next);
        upcase (keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "GKC") {
            iss >> gkc;
        }
        else if (keyword == "GK-RADIUS") {
            getword(record,value,next);
            upcase(value);
            if (value == "VDW") {
               radtyp = "VDW";
            }
            else if (value == "MACROMODEL") {
               radtyp = "MACROMODEL";
            }
            else if (value == "AMOEBA") {
               radtyp = "AMOEBA";
            }
            else if (value == "BONDI") {
               radtyp = "BONDI";
            }
            else if (value == "TOMASI") {
               radtyp = "TOMASI";
            }
            else if (value == "SOLUTE") {
               radtyp = "SOLUTE";
            }
        }
        else if (keyword == "DESCREEN-VDW") {
            getword(record,value,next);
            upcase(value);
            if (value == "TRUE") {
               descreenVDW = true;
            }
        }
        else if (keyword == "DESCREEN-HYDROGEN") {
            getword (record,value,next);
            upcase (value);
            if (value == "FALSE") {
               descreenHydrogen = false;
            }
        }
        else if (keyword == "HCT-SCALE") {
            iss >> dhct;
        }
    }

    // determine the solute atomic radii values to be used
    setrad(radtyp);

    // assign generic value for the HCT overlap scale factor
    for (int i = 0; i < n; i++) {
        shct[i] = dhct;
        rdescr[i] = rsolv[i];
        if (descreenVDW) {
            it = jvdw[i];
            rdescr[i] = 0.5 * radmin[it][it];
        }
        if (!descreenHydrogen) {
            atmnum = atomic[i];
            if (atmnum == 1) shct[i] = 0.;
        }
    }
    if (radtyp == "MACROMODEL") {
        for (int i = 0; i < n; i++) {
            shct[i] = 0.8;
            atmnum = atomic[i];
            if (atmnum == 1) shct[i] = 0.85;
            if (atmnum == 6) shct[i] = 0.72;
            if (atmnum == 7) shct[i] = 0.79;
            if (atmnum == 8) shct[i] = 0.85;
            if (atmnum == 9) shct[i] = 0.88;
            if (atmnum == 15) shct[i] = 0.86;
            if (atmnum == 16) shct[i] = 0.96;
            if (atmnum == 26) shct[i] = 0.88;
        }
    }
}

////////////////////////////////////////////////////
//                                                //
//  kpb  --  assign Poisson-Boltzmann parameters  //
//                                                //
////////////////////////////////////////////////////

// "kpb" assigns parameters needed for the Poisson-Boltzmann
// implicit solvation model implemented via APBS

void kpb()
{

}
//       use sizes
//       use atomid
//       use atoms
//       use bath
//       use couple
//       use gkstuf
//       use inform
//       use iounit
//       use keys
//       use kvdws
//       use math
//       use nonpol
//       use pbstuf
//       use polar
//       use polopt
//       use polpot
//       use potent
//       use ptable
//       use solute
//       implicit none
//       integer i,j
//       integer nx,ny,nz
//       integer maxgrd,next
//       integer pbtyplen,pbsolnlen
//       integer bcfllen,chgmlen
//       integer srfmlen,pbionq
//       integer trimtext
//       real*8 ri,spacing
//       real*8 gx,gy,gz
//       real*8 xcm,ycm,zcm
//       real*8 total,weigh
//       real*8 xmin,xmax,ymin
//       real*8 ymax,zmin,zmax
//       real*8 xlen,ylen,zlen,minlen
//       real*8 pbionc,pbionr
//       character*10 radtyp
//       character*20 keyword
//       character*20 value
//       character*240 record
//       character*240 string
// c
// c
// c     perform dynamic allocation of some global arrays
// c
//       if (allocated(shct))  deallocate (shct)
//       if (allocated(udirs))  deallocate (udirs)
//       if (allocated(udirps))  deallocate (udirps)
//       if (allocated(uinds))  deallocate (uinds)
//       if (allocated(uinps))  deallocate (uinps)
//       if (allocated(uopts))  deallocate (uopts)
//       if (allocated(uoptps))  deallocate (uoptps)
//       allocate (shct(n))
//       allocate (udirs(3,n))
//       allocate (udirps(3,n))
//       allocate (uinds(3,n))
//       allocate (uinps(3,n))
//       if (poltyp == "OPT") {
//         allocate (uopts(0:optorder,3,n))
//         allocate (uoptps(0:optorder,3,n))
//       }
// c
// c     assign some default APBS configuration parameters
// c
//       pbtyp = "LPBE"
//       pbsoln = "MG-MANUAL"
//       radtyp = "SOLUTE"
//       chgm = "SPL4"
//       srfm = "MOL"
//       bcfl = "MDH"
//       kelvin = 298.
//       pdie = 1.
//       sdie = 78.3
//       srad = 0.
//       swin = 0.3
//       sdens = 10.
//       smin = 3.
//       ionn = 0
//       do i = 1, maxion
//          ionc[i] = 0.
//          ionq[i] = 1
//          ionr[i] = 2.
//       }
//       spacing = 0.5
//       maxgrd = 513
// c
// c     compute the position of the center of mass
// c
//       total = 0.
//       xcm = 0.
//       ycm = 0.
//       zcm = 0.
//       do i = 1, n
//          weigh = mass[i]
//          total = total + weigh
//          xcm = xcm + x[i]*weigh
//          ycm = ycm + y[i]*weigh
//          zcm = zcm + z[i]*weigh
//       }
//       xcm = xcm / total
//       ycm = ycm / total
//       zcm = zcm / total
//       gcent(1) = xcm
//       gcent(2) = ycm
//       gcent(3) = zcm
// c
// c     set default APBS grid dimension based on system extent
// c
//       xmin = xcm
//       ymin = ycm
//       zmin = zcm
//       xmax = xcm
//       ymax = ycm
//       zmax = zcm
//       do i = 1, n
//          ri = 1.0
//          xmin = min(xmin,x[i]-ri)
//          ymin = min(ymin,y[i]-ri)
//          zmin = min(zmin,z[i]-ri)
//          xmax = max(xmax,x[i]+ri)
//          ymax = max(ymax,y[i]+ri)
//          zmax = max(zmax,z[i]+ri)
//       }
//       xlen = 2. * (max(xcm-xmin,xmax-xcm)+smin)
//       ylen = 2. * (max(ycm-ymin,ymax-ycm)+smin)
//       zlen = 2. * (max(zcm-zmin,zmax-zcm)+smin)
//       dime(1) = int(xlen/spacing) + 1
//       dime(2) = int(ylen/spacing) + 1
//       dime(3) = int(zlen/spacing) + 1
// c
// c     get any altered APBS parameters from the keyfile
// c
//       do i = 1, nkey
//          next = 1
//          record = keyline[i]
//          call gettext (record,keyword,next)
//          call upcase (keyword)
//          string = record(next:240)
//          if (keyword(1:13) == "APBS-MG-AUTO") {
//             pbsoln = "MG-AUTO"
//          else if (keyword(1:15) == "APBS-MG-MANUAL") {
//             pbsoln = "MG-MANUAL"
//          else if (keyword(1:10) == "APBS-GRID") {
//             nx = dime(1)
//             ny = dime(2)
//             nz = dime(3)
//             read (string,*,err=10,end=10)  nx, ny, nz
//    10       continue
//             if (nx >= 33)  dime(1) = nx
//             if (ny >= 33)  dime(2) = ny
//             if (nz >= 33)  dime(3) = nz
//          else if (keyword(1:11) == "APBS-RADII") {
//             call getword (record,value,next)
//             call upcase (value)
//             if (value(1:3) == "VDW") {
//                radtyp = "VDW"
//             else if (value(1:10) == "MACROMODEL") {
//                radtyp = "MACROMODEL"
//             else if (value(1:5) == "BONDI") {
//                radtyp = "BONDI"
//             else if (value(1:6) == "TOMASI") {
//                radtyp = "TOMASI"
//             else if (value(1:6) == "SOLUTE") {
//                radtyp = "SOLUTE"
//             }
//          else if (keyword(1:11) == "APBS-SDENS") {
//             read (string,*,err=20,end=20)  sdens
//    20       continue
//          else if (keyword(1:10) == "APBS-PDIE") {
//             read (string,*,err=30,end=30)  pdie
//    30       continue
//          else if (keyword(1:10) == "APBS-SDIE") {
//             read (string,*,err=40,end=40)  sdie
//    40       continue
//          else if (keyword(1:10) == "APBS-SRAD") {
//             read (string,*,err=50,end=50)  srad
//    50       continue
//          else if (keyword(1:10) == "APBS-SWIN") {
//             read (string,*,err=60,end=60)  swin
//    60       continue
//          else if (keyword(1:10) == "APBS-SMIN") {
//             read (string,*,err=70,end=70)  smin
//    70       continue
//          else if (keyword(1:7) == "PBTYPE") {
//             call getword (record,value,next)
//             call upcase (value)
//             if (value(1:4) == "LPBE") {
//                pbtyp = "LPBE"
//             else if (value(1:4) == "NPBE") {
//                pbtyp = "NPBE"
//             }
//          else if (keyword(1:10) == "APBS-CHGM") {
//             call getword (record,value,next)
//             call upcase (value)
//             if (value(1:4) == "SPL0") {
//                chgm = "SPL0"
//             else if (value(1:4) == "SPL2") {
//                chgm = "SPL2"
//             else if (value(1:4) == "SPL4") {
//                chgm = "SPL4"
//             }
//          else if (keyword(1:10) == "APBS-SRFM") {
//             call getword (record,value,next)
//             call upcase (value)
//             if (value(1:3) == "MOL") {
//                srfm = "MOL"
//             else if (value(1:4) == "SMOL") {
//                srfm = "SMOL"
//             else if (value(1:4) == "SPL2") {
//                srfm = "SPL2"
//             else if (value(1:4) == "SPL4") {
//                srfm = "SPL4"
//             }
//          else if (keyword(1:10) == "APBS-BCFL") {
//             call getword (record,value,next)
//             call upcase (value)
//             if (value(1:3) == "ZERO") {
//                bcfl = "ZERO"
//             else if (value(1:3) == "MDH") {
//                bcfl = "MDH"
//             else if (value(1:3) == "SDH") {
//                bcfl = "SDH"
//             }
//          else if (keyword(1:9) == "APBS-ION") {
//             pbionc = 0.
//             pbionq = 1
//             pbionr = 2.
//             read (string,*,err=80,end=80)  pbionq,pbionc,pbionr
//    80       continue
//             if (pbionq!=0 and pbionc>=0.
//      &             and pbionr>=0.) {
//                ionn = ionn + 1
//                ionc(ionn) = pbionc
//                ionq(ionn) = pbionq
//                ionr(ionn) = pbionr
//             }
//          }
//       }
// c
// c     set APBS grid spacing for the chosen grid dimension
// c
//       xlen = 2. * (max(xcm-xmin,xmax-xcm)+smin)
//       ylen = 2. * (max(ycm-ymin,ymax-ycm)+smin)
//       zlen = 2. * (max(zcm-zmin,zmax-zcm)+smin)
//       grid(1) = xlen / dime(1)
//       grid(2) = ylen / dime(2)
//       grid(3) = zlen / dime(3)
// c
// c     grid spacing must be equal to maintain traceless quadrupoles
// c
//       grid(1) = min(grid(1),grid(2),grid(3))
//       grid(2) = grid(1)
//       grid(3) = grid(1)
// c
// c     set the grid dimensions to the smallest multiples of 32
// c
//       dime(1) = 33
//       dime(2) = 33
//       dime(3) = 33
// c
// c     use minimum side length to maintain equal grid spacing
// c
//       minlen = min(xlen,ylen,zlen)
//       do while (grid(1)*dime(1) < minlen)
//          dime(1) = dime(1) + 32
//       }
//       do while (grid(2)*dime(2) < minlen)
//          dime(2) = dime(2) + 32
//       }
//       do while (grid(3)*dime(3) < minlen)
//          dime(3) = dime(3) + 32
//       }
// c
// c     limit the grid dimensions and recompute the grid spacing
// c
//       dime(1) = min(dime(1),maxgrd)
//       dime(2) = min(dime(2),maxgrd)
//       dime(3) = min(dime(3),maxgrd)
//       grid(1) = xlen / dime(1)
//       grid(2) = ylen / dime(2)
//       grid(3) = zlen / dime(3)
// c
// c     grid spacing must be equal to maintain traceless quadrupoles
// c
//       grid(1) = max(grid(1),grid(2),grid(3))
//       grid(2) = grid(1)
//       grid(3) = grid(1)
// c
// c     if this is an "mg-auto" (focusing) calculation, set the
// c     fine grid to the default size, and the coarse grid to
// c     twice its original size; currently, all energies and
// c     forces need to be evaluated at the same resolution
// c
//       if (pbsoln == "MG-AUTO") {
//          fgrid(1) = grid(1)
//          fgrid(2) = grid(2)
//          fgrid(3) = grid(3)
//          fgcent(1) = gcent(1)
//          fgcent(2) = gcent(2)
//          fgcent(3) = gcent(3)
//          cgrid(1) = 2. * grid(1)
//          cgrid(2) = 2. * grid(2)
//          cgrid(3) = 2. * grid(3)
//       }
// c
// c     get any custom APBS grid parameters from the keyfile
// c
//       do i = 1, nkey
//          next = 1
//          record = keyline[i]
//          call gettext (record,keyword,next)
//          call upcase (keyword)
//          string = record(next:240)
//          if (keyword(1:10) == "APBS-DIME") {
//             read (string,*,err=90,end=90)  nx,ny,nz
//             dime(1) = nx
//             dime(2) = ny
//             dime(3) = nz
//    90       continue
//             do j = 1, 3
//                if (mod(dime(j),32) != 1) {
//                   dime(j) = 32*(1+(dime(j)-1)/32) + 1
//                }
//             }
//          else if (keyword(1:11) == "APBS-AGRID") {
//             read (string,*,err=100,end=100)  gx,gy,gz
//             grid(1) = gx
//             grid(2) = gy
//             grid(3) = gz
//   100       continue
//          else if (keyword(1:11) == "APBS-CGRID") {
//             read (string,*,err=110,end=110)  gx,gy,gz
//             cgrid(1) = gx
//             cgrid(2) = gy
//             cgrid(3) = gz
//   110       continue
//          else if (keyword(1:11) == "APBS-FGRID") {
//             read (string,*,err=120,end=120)  gx,gy,gz
//             fgrid(1) = gx
//             fgrid(2) = gy
//             fgrid(3) = gz
//   120       continue
//          else if (keyword(1:11) == "APBS-GCENT") {
//             read (string,*,err=130,end=130)  gx,gy,gz
//             gcent(1) = gx
//             gcent(2) = gy
//             gcent(3) = gz
//   130       continue
//          else if (keyword(1:12) == "APBS-CGCENT") {
//             read (string,*,err=140,end=140)  gx,gy,gz
//             cgcent(1) = gx
//             cgcent(2) = gy
//             cgcent(3) = gz
//   140       continue
//          else if (keyword(1:12) == "APBS-FGCENT") {
//             read (string,*,err=150,end=150)  gx,gy,gz
//             fgcent(1) = gx
//             fgcent(2) = gy
//             fgcent(3) = gz
//   150       continue
//          }
//       }
// c
// c     determine the solute atomic radii values to be used
// c
//       call setrad (radtyp)
// c
// c     assign generic value for the HCT overlap scale factor
// c
//       do i = 1, n
//          shct[i] = 0.69
//       }
// c
// c     determine the length of the character arguments
// c
//       pbtyplen = trimtext (pbtyp)
//       pbsolnlen = trimtext (pbsoln)
//       bcfllen = trimtext (bcfl)
//       chgmlen = trimtext (chgm)
//       srfmlen = trimtext (srfm)
// c
// c     make call needed to initialize the APBS calculation
// c
//       call apbsinitial (dime,grid,gcent,cgrid,cgcent,fgrid,fgcent,
//      &                  pdie,sdie,srad,swin,sdens,kelvin,ionn,ionc,
//      &                  ionq,ionr,pbtyp,pbtyplen,pbsoln,pbsolnlen,
//      &                  bcfl,bcfllen,chgm,chgmlen,srfm,srfmlen)
// c
// c     print out the APBS grid dimensions and spacing
// c
//       if (verbose) {
//          write (iout,160)  (dime[i],i=1,3),grid(1)
//   160    format (/," APBS Grid Dimensions and Spacing :",
//      &           //,10x,3i8,10x,f10.4)
//       }

////////////////////////////////////////////////////
//                                                //
//  knp  --  assign cavity-dispersion parameters  //
//                                                //
////////////////////////////////////////////////////

// "knp" initializes parameters needed for the cavity-plus-
// dispersion nonpolar implicit solvation model

void knp()
{
    int i,next;
    double cross,ah,ao;
    double rmini,epsi;
    double rmixh,rmixh3;
    double rmixh7,emixh;
    double rmixo,rmixo3;
    double rmixo7,emixo;
    double ri,ri3,ri7,ri11;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // set default values for solvent pressure and surface tension
    solvprs = 0.0334;
    surften = 0.103;

    // get any altered surface tension value from keyfile
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "SOLVENT-PRESSURE") {
            iss >> solvprs;
        }
        else if (keyword == "SURFACE-TENSION") {
            iss >> surften;
        }
    }

    // set switching function values for pressure and tension
    // cross = 9.251 = 3.0 * 0.103 / 0.0334
    cross = 3. * surften / solvprs;
    spcut = cross - 3.5;
    spoff = cross + 3.5;

    // The SASA term is switched on 0.2 Angtroms after 
    // the cross-over point to give a smooth transition
    stcut = cross + 3.5 + 0.2;
    stoff = cross - 3.5 + 0.2;

    // perform dynamic allocation of some global arrays
    if (asolv.size() != 0) asolv.resize(0);
    if (radcav.size() != 0) radcav.resize(0);
    if (raddsp.size() != 0) raddsp.resize(0);
    if (epsdsp.size() != 0) epsdsp.resize(0);
    if (cdsp.size() != 0) cdsp.resize(0);
    asolv.resize(n);
    radcav.resize(n);
    raddsp.resize(n);
    epsdsp.resize(n);
    cdsp.resize(n);

    // assign surface area factors for nonpolar solvation
    for (int i = 0; i < n; i++) {
        asolv[i] = surften;
    }

    // set cavity and dispersion radii for nonpolar solvation
    for (int i = 0; i < n; i++) {
        if (vdwindex == "CLASS") {
            radcav[i] = rad[atomClass[i]] + cavoff;
            raddsp[i] = rad[atomClass[i]];
            epsdsp[i] = eps[atomClass[i]];
        }
        else {
            radcav[i] = rad[type[i]] + cavoff;
            raddsp[i] = rad[type[i]];
            epsdsp[i] = eps[type[i]];
        }
    }

    // compute maximum dispersion energies for each atom
    for (int i = 0; i < n; i++) {
        epsi = epsdsp[i];
        rmini = raddsp[i];
        if (rmini>0. and epsi>0.) {
            emixo = 4. * epso * epsi / std::pow((std::sqrt(epso)+std::sqrt(epsi)),2);
            rmixo = 2. * (std::pow(rmino,3)+std::pow(rmini,3)) / (std::pow(rmino,2)+std::pow(rmini,2));
            rmixo3 = std::pow(rmixo,3);
            rmixo7 = std::pow(rmixo,7);
            ao = emixo * rmixo7;
            emixh = 4. * epsh * epsi / std::pow((std::sqrt(epsh)+std::sqrt(epsi)),2);
            rmixh = 2. * (std::pow(rminh,3)+std::pow(rmini,3)) / (std::pow(rminh,2)+std::pow(rmini,2));
            rmixh3 = std::pow(rmixh,3);
            rmixh7 = std::pow(rmixh,7);
            ah = emixh * rmixh7;
            ri = 0.5*rmixh + dspoff;
            ri3 = std::pow(ri,3);
            ri7 = std::pow(ri,7);
            ri11 = std::pow(ri,11);
            if (ri < rmixh) {
                cdsp[i] = -4.*pi*emixh*(rmixh3-ri3)/3.;
                cdsp[i] = cdsp[i] - emixh*18./11.*rmixh3*pi;
            }
            else {
                cdsp[i] = 2.*pi*(2.*rmixh7-11.*ri7)*ah;
                cdsp[i] = cdsp[i] / (11.*ri11);
            }
            cdsp[i] = 2. * cdsp[i];
            ri = 0.5*rmixo + dspoff;
            ri3 = std::pow(ri,3);
            ri7 = std::pow(ri,7);
            ri11 = std::pow(ri,11);
            if (ri < rmixo) {
                cdsp[i] = cdsp[i] - 4.*pi*emixo*(rmixo3-ri3)/3.;
                cdsp[i] = cdsp[i] - emixo*18./11.*rmixo3*pi;
            }
            else {
                cdsp[i] = cdsp[i] + 2.*pi*(2.*rmixo7-11.*ri7) * ao/(11.*ri11);
            }
        }
        cdsp[i] = slevy * awater * cdsp[i];
    }
}

////////////////////////////////////////////////////
//                                                //
//  khpmf  --  assign hydrophobic PMF parameters  //
//                                                //
////////////////////////////////////////////////////

//  "khpmf" initializes parameters needed for the hydrophobic
//  potential of mean force nonpolar implicit solvation model
// 
//  literature reference:
// 
//  M. S. Lin, N. L. Fawzi and T. Head-Gordon, "Hydrophobic
//  Potential of Mean Force as a Solvation Function for Protein
//  Structure Prediction", Structure, 15, 727-740 (2007)

void khpmf()
{

}
//       use sizes
//       use atomid
//       use atoms
//       use couple
//       use hpmf
//       use ptable
//       implicit none
//       integer i,j,k
//       integer nh,atn
//       logical keep
// c
// c
// c     perform dynamic allocation of some global arrays
// c
//       if (allocated(ipmf))  deallocate (ipmf)
//       if (allocated(rpmf))  deallocate (rpmf)
//       if (allocated(acsa))  deallocate (acsa)
//       allocate (ipmf(n))
//       allocate (rpmf(n))
//       allocate (acsa(n))
// c
// c     get carbons for PMF and set surface area screening values
// c
//       npmf = 0
//       do i = 1, n
//          if (atomic[i] == 6) {
//             keep = true
//             nh = 0
//             if (n12[i] <= 2)  keep = false
//             do j = 1, n12[i]
//                k = i12(j,i)
//                if (atomic[k] == 1)  nh = nh + 1
//                if (n12[i]==3 and atomic[k]==8)  keep = false
//             }
//             if (keep) {
//                npmf = npmf + 1
//                ipmf(npmf) = i
//                acsa[i] = 1.
//                if (n12[i]==3 and nh==0)  acsa[i] = 1.554
//                if (n12[i]==3 and nh==1)  acsa[i] = 1.073
//                if (n12[i]==4 and nh==1)  acsa[i] = 1.276
//                if (n12[i]==4 and nh==2)  acsa[i] = 1.045
//                if (n12[i]==4 and nh==3)  acsa[i] = 0.88
//                acsa[i] = acsa[i] * safact/acsurf
//             }
//          }
//       }
// c
// c     assign HPMF atomic radii from consensus vdw values
// c
//       do i = 1, n
//          rpmf[i] = 1.
//          atn = atomic[i]
//          if (atn == 0) {
//             rpmf[i] = 0.0
//          else
//             rpmf[i] = vdwrad(atn)
//          }
//          if (atn == 5)  rpmf[i] = 1.8
//          if (atn == 8)  rpmf[i] = 1.5
//          if (atn == 35)  rpmf[i] = 1.85
//       }

/////////////////////////////////////////////////////
//                                                 //
//  setrad  --  assign solute radii for PB and GK  //
//                                                 //
/////////////////////////////////////////////////////

// "setrad" chooses a set of atomic radii to solute atoms for use
// during Poission-Boltzmann and Generalized Kirkwood implicit
// solvation calculations

void setrad(std::string radtyp)
{
    int i,j,k,l,m;
    int atmnum;
    double rscale;
    double offset;

    // assign default solute radii from consensus vdw values
    for (int i = 0; i < n; i++) {
        atmnum = atomic[i];
        if (atmnum == 0)  rsolv[i] = 0.;
        rsolv[i] = vdwrad[atmnum-1];
    }

    // assign solute atomic radii from force field vdw values
    if (radtyp == "VDW") {
        for (int i = 0; i < n; i++) {
            k = jvdw[i];
            rsolv[i] = 2.;
            if (k != -1) {
                rsolv[i] = 0.5 * radmin[k][k];
            }
        }
    }    

    // assign solute radii from parametrized solvation values
    else if (radtyp == "SOLUTE") {
        if (solvtyp == "GK") {     
            for (int i = 0; i < n; i++) {
                if (type[i] != -1) {
                    if (gkr[type[i]] != 0.) {
                        rsolv[i] = gkr[type[i]];
                    }
                }
            }
        }
        else if (solvtyp == "PB") {
            for (int i = 0; i < n; i++) {
                if (type[i] != -1) {
                    if (pbr[type[i]] != 0.) {
                        rsolv[i] = pbr[type[i]];
                    }
                }
            }
        }
    }

    // assign solute atomic radii adapted from Macromodel
    else if (radtyp == "MACROMODEL") {
        for (int i = 0; i < n; i++) {
            atmnum = atomic[i];
            if (atmnum == 0) rsolv[i] = 0.;
            rsolv[i] = vdwrad[atmnum-1];
            if (atmnum == 1) {
                rsolv[i] = 1.25;
                k = i12[i][0];
                if (atomic[k] == 7) rsolv[i] = 1.15;
                if (atomic[k] == 8) rsolv[i] = 1.05;
            }
            else if (atmnum == 3) {
                rsolv[i] = 1.432;
            }
            else if (atmnum == 6) {
                rsolv[i] = 1.9;
                if (n12[i] == 3) rsolv[i] = 1.875;
                if (n12[i] == 2) rsolv[i] = 1.825;
            }
            else if (atmnum == 7) {
                rsolv[i] = 1.7063;
                if (n12[i] == 4) rsolv[i] = 1.625;
                if (n12[i] == 1) rsolv[i] = 1.6;
            }
            else if (atmnum == 8) {
                rsolv[i] = 1.535;
                if (n12[i] == 1) rsolv[i] = 1.48;
            }
            else if (atmnum == 9) {
                rsolv[i] = 1.47;
            }
            else if (atmnum == 10) {
                rsolv[i] = 1.39;
            }
            else if (atmnum == 11) {
                rsolv[i] = 1.992;
            }
            else if (atmnum == 12) {
                rsolv[i] = 1.7;
            }
            else if (atmnum == 14) {
                rsolv[i] = 1.8;
            }
            else if (atmnum == 15) {
                rsolv[i] = 1.87;
            }
            else if (atmnum == 16) {
                rsolv[i] = 1.775;
            }
            else if (atmnum == 17) {
                rsolv[i] = 1.735;
            }
            else if (atmnum == 18) {
                rsolv[i] = 1.7;
            }
            else if (atmnum == 19) {
                rsolv[i] = 2.123;
            }
            else if (atmnum == 20) {
                rsolv[i] = 1.817;
            }
            else if (atmnum == 35) {
                rsolv[i] = 1.9;
            }
            else if (atmnum == 36) {
                rsolv[i] = 1.812;
            }
            else if (atmnum == 37) {
                rsolv[i] = 2.26;
            }
            else if (atmnum == 53) {
                rsolv[i] = 2.1;
            }
            else if (atmnum == 54) {
                rsolv[i] = 1.967;
            }
            else if (atmnum == 55) {
                rsolv[i] = 2.507;
            }
            else if (atmnum == 56) {
                rsolv[i] = 2.188;
            }
        }
    }

    // assign solute atomic radii as modified Bondi values
    else if (radtyp == "AMOEBA") {
        for (int i = 0; i < n; i++) {
            atmnum = atomic[i];
            if (atmnum == 0) rsolv[i] = 0.;
            rsolv[i] = vdwrad[atmnum-1];
            if (atmnum == 1) {
                rsolv[i] = 1.32;
                k = i12[i][0];
                if (atomic[k] == 7) rsolv[i] = 1.1;
                if (atomic[k] == 8) rsolv[i] = 1.05;
            }
            if (atmnum == 3) rsolv[i] = 1.5;
            if (atmnum == 6) {
                rsolv[i] = 2.0;
                if (n12[i] == 3) rsolv[i] = 2.05;
                if (n12[i] == 4) {
                    for (int j = 0; j < n12[i]; j++) {
                        k = i12[i][j];
                        if (atomic[k] == 7) rsolv[i] = 1.75;
                        if (atomic[k] == 8) rsolv[i] = 1.75;
                    }
                }
            }
            if (atmnum == 7) {
                rsolv[i] = 1.6;
            }
            if (atmnum == 8) {
                rsolv[i] = 1.55;
                if (n12[i] == 2) rsolv[i] = 1.45;
            }
        }
    }

    // make Tomasi-style modifications to the solute radii values
    else if (radtyp == "TOMASI") {
        for (int i = 0; i < n; i++) {
            offset = 0.;
            atmnum = atomic[i];
            if (atomic[i] == 1) {
                for (int j = 0; j < n12[i]; j++) {
                    k = i12[i][j];
                    if (atomic[k] == 6) {
                        for (int l = 0; l < n12[k]; l++) {
                            m = i12[k][l];
                            if (atomic[m] == 7) offset = -0.05;
                            if (atomic[m] == 8) offset = -0.1;
                        }
                    }
                    if (atomic[k] == 7) offset = -0.25;
                    if (atomic[k] == 8) offset = -0.4;
                    if (atomic[k] == 16) offset = -0.1;
                }
            }
            else if (atomic[i] == 6) {
                if (n12[i] == 4) offset = 0.05;
                if (n12[i] == 3) offset = 0.02;
                if (n12[i] == 2) offset = -0.03;
                for (int j = 0; j < n12[i]; j++) {
                    k = i12[i][j];
                    if (atomic[k] == 6) offset = offset - 0.07;
                }
                for (int j = 0; j < n12[i]; j++) {
                    k = i12[i][j];
                    if (atomic[k]==7 and n12[k]==4) offset = -0.2;
                    if (atomic[k]==7 and n12[k]==3) offset = -0.25;
                    if (atomic[k] == 8) offset = -0.2;
                }
            }
            else if (atomic[i] == 7) {
                if (n12[i] == 3) {
                    offset = -0.1;
                    for (int j = 0; j < n12[i]; j++) {
                        k = i12[i][j];
                        if (atomic[k] == 6) offset = offset - 0.24;
                    }
                }
                else {
                    offset = -0.2;
                    for (int j = 0; j < n12[i]; j++) {
                        k = i12[i][j];
                        if (atomic[k] == 6) offset = offset - 0.16;
                    }
                }
            }
            else if (atomic[i] == 8) {
                if (n12[i] == 2) {
                    offset = -0.21;
                    for (int j = 0; j < n12[i]; j++) {
                        k = i12[i][j];
                        if (atomic[k] == 6) offset = -0.36;
                    }
                }
                else {
                    offset = -0.25;
                }
            }
            else if (atomic[i] == 16) {
                offset = -0.03;
                for (int j = 0; j < n12[i]; j++) {
                    k = i12[i][j];
                    if (atomic[k] == 6) offset = offset - 0.1;
                }
            }
            rsolv[i] = rsolv[i] + offset;
        }
    }

    // apply an overall scale factor to the solute atomic radii
    rscale = 1.;
    if (radtyp == "MACROMODEL") rscale = 1.15;
    if (radtyp == "BONDI") rscale = 1.21;
    if (radtyp == "TOMASI") rscale = 1.47;
    for (int i = 0; i < n; i++) {
        rsolv[i] = rsolv[i] * rscale;
    }
}
}
