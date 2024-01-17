// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "chgpen.h"
#include "chgtrn.h"
#include "chkpole.h"
#include "expol.h"
#include "inform.h"
#include "kchgtrn.h"
#include "kctrn.h"
#include "keys.h"
#include "mplpot.h"
#include "mpole.h"
#include "polar.h"
#include "polpot.h"
#include "potent.h"
#include "sizes.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  kchgtrn  --  charge transfer term assignment  //
//                                                //
////////////////////////////////////////////////////

// "kchgtrn" assigns charge magnitude and damping parameters for
// charge transfer interactions and processes any new or changed
// values for these parameters

// TODO
void kchgtrn()
{
    int i,k;
    int ia,ic,next;
    real chtrn,actrn;
    bool header;
    std::string keyword;
    std::string record;
    std::string string;

// c     process keywords containing charge transfer parameters
// c
//       header = .true.
//       do i = 1, nkey
//          next = 1
//          record = keyline(i)
//          call gettext (record,keyword,next)
//          call upcase (keyword)
//          if (keyword(1:7) .eq. 'CHGTRN ') then
//             k = 0
//             chtrn = 0.0d0
//             actrn = 0.0d0
//             call getnumb (record,k,next)
//             string = record(next:240)
//             read (string,*,err=10,end=10) chtrn,actrn
//    10       continue
//             if (k .gt. 0) then
//                if (header .and. .not.silent) then
//                   header = .false.
//                   write (iout,20)
//    20             format (/,' Additional Charge Transfer',
//      &                       ' Parameters :',
//      &                    //,5x,'Atom Class',13x,'Charge',11x,'Damp',/)
//                end if
//                if (k .le. maxclass) then
//                   ctchg(k) = chtrn
//                   ctdmp(k) = actrn
//                   if (.not. silent) then
//                      write (iout,30) k,chtrn,actrn
//    30                format (6x,i6,7x,f15.4,f15.4)
//                   end if
//                else
//                   write (iout,40)
//    40             format (/,' KCHGTRN  --  Too many Charge',
//      &                       ' Transfer Parameters')
//                   abort = .true.
//                end if
//             end if
//          end if
//       end do

    // perform dynamic allocation of some global arrays
    if (chgct.size() != 0) chgct.resize(0);
    if (dmpct.size() != 0) dmpct.resize(0);
    chgct.resize(n);
    dmpct.resize(n);
// c
// c     assign the charge transfer charge and alpha parameters 
// c
//       nct = n
//       do i = 1, n
//          ic = class(i)
//          chgct(i) = ctchg(ic)
//          dmpct(i) = ctdmp(ic)
//       end do
// c
// c     process keywords containing atom specific charge transfer
// c
//       header = .true.
//       do i = 1, nkey
//          next = 1
//          record = keyline(i)
//          call gettext (record,keyword,next)
//          call upcase (keyword)
//          if (keyword(1:7) .eq. 'CHGTRN ') then
//             ia = 0
//             chtrn = 0.0d0
//             actrn = 0.0d0
//             string = record(next:240)
//             read (string,*,err=70,end=70) ia,chtrn,actrn
//             if (ia.lt.0 .and. ia.ge.-n) then
//                ia = -ia
//                if (header .and. .not.silent) then
//                   header = .false.
//                   write (iout,50)
//    50             format (/,' Additional Charge Transfer Values',
//      &                       ' for Specific Atoms :',
//      &                    //,8x,'Atom',16x,'Charge',11x,'Damp',/)
//                end if
//                if (.not. silent) then
//                   write (iout,60) ia,chtrn,actrn
//    60             format (6x,i6,7x,f15.4,f15.4)
//                end if
//                chgct(ia) = chtrn
//                dmpct(ia) = actrn
//             end if
//    70       continue
//          end if
//       end do

    // remove zero or undefined electrostatic sites from the list
    if (use_chgtrn) {
        npole = 0;
        ncp = 0;
        npolar = 0;
        nexpol = 0;
        nct = 0;
        for (int i = 0; i < n; i++) {
            if (polarity[i] == 0.) douind[i] = false;
            if (polsiz[i]!=0 or polarity[i]!=0. or chgct[i]!=0. or dmpct[i]!=0.) {
                ipole[npole] = i;
                pollist[i] = npole;
                mono0[i] = pole[i][0];
                if (palpha[i] != 0.) ncp++;
                if (polarity[i] != 0.) {
                    ipolar[npolar] = npole;
                    npolar++;
                    douind[i] = true;
                }
                if (tholed[i] != 0.) use_tholed = true;
                // if (kpep[i] != 0.) nexpol++;
                if (chgct[i]!=0. or dmpct[i]!=0.) nct++;
                npole++;
            }
        }
    }

    // test multipoles at chiral sites and invert if necessary
    if (use_chgtrn) chkpole();

    // turn off individual electrostatic potentials if not used
    if (npole == 0) use_mpole = false;
    if (npolar == 0) use_polar = false;
    if (ncp != 0) use_chgpen = true;
    if (ncp != 0) use_thole = false;
    if (use_tholed) use_thole = true;
    if (nexpol != 0) use_expol = true;
    if (nct == 0) use_chgtrn = false;
}
}
