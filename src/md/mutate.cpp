// Author: Moses KJ Chung
// Year:   2024

#include "atomid.h"
#include "atoms.h"
#include "bndstr.h"
#include "gettext.h"
#include "inform.h"
#include "katoms.h"
#include "keys.h"
#include "mutant.h"
#include "mutate.h"
#include "potent.h"
#include "upcase.h"
#include <sstream>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  mutate  --  set parameters for hybrid system  //
//                                                //
////////////////////////////////////////////////////

// "mutate" constructs the hybrid hamiltonian for a specified
// initial state, final state and mutation parameter "lambda"
//
// note torsional and most electrostatics terms apply "lambda"
// by directly scaling parameters, while vdw and repulsion energy
// terms use soft core functions from the references cited below
//
// literature references:
//
// T. Steinbrecher, D. L. Mobley and D. A. Case, "Nonlinear Scaling
// Schemes for Lennard-Jones Interactions in Free Energy
// Calculations", Journal of Chemical Physics, 127, 214108 (2007)
//
// D. Jiao, P. A. Golubkov, T. A. Darden and P. Ren, "Calculation
// of Protein-Ligand Binding Free Energy by Using a Polarizable
// Potential", PNAS, 105, 6290-6295 (2008)

void mutate()
{
    int j,k,ihyb;
    int it0,it1;
    int next,size;
    int ntbnd;
    std::vector<int> list;
    std::vector<std::vector<int>> itbnd;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // perform dynamic allocation of some global arrays
    imut.allocate(n);
    type0.allocate(n);
    class0.allocate(n);
    type1.allocate(n);
    class1.allocate(n);
    mut.allocate(n);

    // perform dynamic allocation of some local arrays
    size = 40;
    list.resize(size);
    itbnd.resize(nbond,std::vector<int>(2));

    // set defaults for lambda perturbation scaling values
    lambda = 1.;
    vlambda = 1.;
    elambda = 1.;
    tlambda = 1.;

    // set defaults for vdw coupling type and soft core vdw
    vcouple = 0;
    scexp = 5.;
    scalpha = 0.7;

    // zero out number of hybrid atoms and mutated torsions
    nmut = 0;
    for (int i = 0; i < n; i++) {
        mut[i] = false;
    }
    ntbnd = 0;
    for (int i = 0; i < nbond; i++) {
        itbnd[i][0] = -1;
        itbnd[i][1] = -1;
    }

    // search keywords for free energy perturbation options
    for (int i = 0; i < nkey; i++) {
        next = 0;
        record = keyline[i];
        gettext(record,keyword,next);
        upcase(keyword);
        if (keyword == "LAMBDA") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> lambda;
        }
        else if (keyword == "VDW-LAMBDA") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> vlambda;
        }
        else if (keyword == "ELE-LAMBDA") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> elambda;
        }
        else if (keyword == "TORS-LAMBDA") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> tlambda;
        }
        else if (keyword == "VDW-ANNIHILATE") {
            vcouple = 1;
        }
        else if (keyword == "MUTATE") {
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            if (!(iss >> ihyb >> it0 >> it1)) break;
            ihyb--;
            it0--;
            it1--;
            imut[nmut] = ihyb;
            mut[ihyb] = true;
            type0[nmut] = it0;
            type1[nmut] = it1;
            class0[nmut] = atmcls[it0];
            class1[nmut] = atmcls[it1];
            nmut++;
        }
        else if (keyword == "LIGAND") {
            for (int k = 0; k < size; k++) {
                list[k] = 0;
            }
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            k = 0;
            while (iss >> list[k]) k++;
            k = 0;
            while (list[k] != 0) {
                if (list[k] > 0) {
                    j = list[k]-1;
                    imut[nmut] = j;
                    mut[j] = true;
                    type0[nmut] = -1;
                    type1[nmut] = type[j];
                    class0[nmut] = -1;
                    class1[nmut] = atomClass[j];
                    nmut++;
                    k++;
                }
                else {
                    int start = std::abs(list[k])-1;
                    int end = std::abs(list[k+1]);
                    for (int j = start; j < end; j++) {
                        imut[nmut] = j;
                        mut[j] = true;
                        type0[nmut] = -1;
                        type1[nmut] = type[i];
                        class0[nmut] = -1;
                        class1[nmut] = atomClass[i];
                        nmut++;
                    }
                    k += 2;
                }
            }
        }
//         else if (keyword == "ROTATABLE-BOND") {
//             do k = 1, size
//                 list(k) = 0
//             end do
//             string = record(next:240)
//             read (string,*,err=20,end=20)  (list(k),k=1,size)
//    20       continue
//             k = 1
//             do while (list(k) != 0)
//                 ntbnd = ntbnd + 1
//                 itbnd(1,ntbnd) = list(k)
//                 itbnd(2,ntbnd) = list(k+1)
//                 k = k + 2
//             end do
//         }
    }
}
}

// c

// c
// c     scale electrostatic parameter values based on lambda
// c
//       if (elambda>=0.0d0 and elambda<1.0d0) then
//          call altelec
//       end if
// c
// c     scale torsional parameter values based on lambda
// c
//       if (tlambda>=0.0d0 and tlambda<1.0d0) then
//          if (ntbnd != 0)  call alttors (ntbnd,itbnd)
//       end if
// c
// c     turn off hybrid potentials if no sites are mutated
// c
//       use_mutate = .true.
//       if (nmut == 0)  use_mutate = .false.
// c
// c     write status of current hybrid potential lambda values
// c
//       if (use_mutate and !silent) then
//          write (iout,40)  vlambda
//    40    format (/,' Free Energy Perturbation :',f15.3,
//      &              ' Lambda for van der Waals')
//          write (iout,50)  elambda
//    50    format (' Free Energy Perturbation :',f15.3,
//      &              ' Lambda for Electrostatics')
//          write (iout,60)  tlambda
//    60    format (' Free Energy Perturbation :',f15.3,
//      &              ' Lambda for Torsional Angles')
//       end if
// c
// c     perform deallocation of some local arrays
// c
//       deallocate (list)
//       deallocate (itbnd)
//       return
//       end
// c
// c
// c     ################################################################
// c     ##                                                            ##
// c     ##  subroutine altelec  --  mutated electrostatic parameters  ##
// c     ##                                                            ##
// c     ################################################################
// c
// c
// c     "altelec" constructs mutated electrostatic parameters based
// c     on the lambda mutation parameter "elambda"
// c
// c     note charge transfer electrostatics is not treated by parameter
// c     scaling due to the functional form used, and must be done via
// c     modification of pairwise energy terms in the potential routines
// c
// c
//       subroutine altelec
//       use angbnd
//       use atoms
//       use bndstr
//       use cflux
//       use charge
//       use chgpen
//       use dipole
//       use mplpot
//       use mpole
//       use mutant
//       use polar
//       use potent
//       implicit none
//       integer i,j,k
//       integer k1,k2
//       integer ia,ib,ic
// c
// c
// c     set scaled parameters for partial charge models
// c
//       if (use_charge) then
//          do i = 1, nion
//             k = iion(i)
//             if (mut(k)) then
//                pchg(k) = pchg(k) * elambda
//             end if
//             pchg0(k) = pchg(k)
//          end do
//       end if
// c
// c     set scaled parameters for bond dipole models
// c
//       if (use_dipole) then
//          do i = 1, ndipole
//             k1 = idpl(1,i)
//             k2 = idpl(2,i)
//             if (mut(k1) or mut(k2)) then
//                bdpl(i) = bdpl(i) * elambda
//             end if
//          end do
//       end if
// c
// c     set scaled parameters for atomic multipole models
// c
//       if (use_mpole) then
//          do i = 1, npole
//             k = ipole(i)
//             if (mut(k)) then
//                do j = 1, 13
//                   pole(j,k) = pole(j,k) * elambda
//                end do
//                mono0(k) = pole(1,k)
//                if (use_chgpen) then
//                   pcore(k) = pcore(k) * elambda
//                   pval(k) = pval(k) * elambda
//                   pval0(k) = pval(k)
//                end if
//             end if
//          end do
//       end if
// c
// c     set scaled parameters for atomic polarizability models
// c
//       if (use_polar) then
//          do i = 1, npole
//             k = ipole(i)
//             if (mut(k)) then
//                polarity(k) = polarity(k) * elambda
//                if (elambda == 0.0d0)  douind(k) = .false.
//             end if
//          end do
//       end if
// c
// c     set scaled parameters for bond stretch charge flux
// c
//       if (use_chgflx) then
//          do i = 1, nbond
//             ia = ibnd(1,i)
//             ib = ibnd(2,i)
//             if (mut(ia) and mut(ib)) then
//                bflx(i) = bflx(i) * elambda
//             end if
//          end do
//       end if
// c
// c     set scaled parameters for angle bend charge flux
// c
//       if (use_chgflx) then
//          do i = 1, nangle
//             ia = iang(1,i)
//             ib = iang(2,i)
//             ic = iang(3,i)
//             if (mut(ia) and mut(ib) and mut(ic)) then
//                aflx(1,i) = aflx(1,i) * elambda
//                aflx(2,i) = aflx(2,i) * elambda
//                abflx(1,i) = abflx(1,i) * elambda
//                abflx(2,i) = abflx(2,i) * elambda
//             end if
//          end do
//       end if
//       return
//       end
// c
// c
// c     ############################################################
// c     ##                                                        ##
// c     ##  subroutine alttors  --  mutated torsional parameters  ##
// c     ##                                                        ##
// c     ############################################################
// c
// c
// c     "alttors" constructs mutated torsional parameters based
// c     on the lambda mutation parameter "tlambda"
// c
// c
//       subroutine alttors (ntbnd,itbnd)
//       use mutant
//       use potent
//       use tors
//       implicit none
//       integer i,j
//       integer ia,ib,ic,id
//       integer kb,kc
//       integer ntbnd
//       integer itbnd(2,*)
// c
// c
// c     set scaled parameters for specified rotatable bonds
// c
//       if (use_tors) then
//          do i = 1, ntors
//             ia = itors(1,i)
//             ib = itors(2,i)
//             ic = itors(3,i)
//             id = itors(4,i)
//             if (mut(ia) and mut(ib) and mut(ic) and mut(id)) then
//                do j = 1, ntbnd
//                   kb = itbnd(1,j)
//                   kc = itbnd(2,j)
//                   if ((kb==ib and kc==ic) or
//      &                (kb==ic and kc==ib)) then
//                      tors1(1,i) = tors1(1,i) * tlambda
//                      tors2(1,i) = tors2(1,i) * tlambda
//                      tors3(1,i) = tors3(1,i) * tlambda
//                      tors4(1,i) = tors4(1,i) * tlambda
//                      tors5(1,i) = tors5(1,i) * tlambda
//                      tors6(1,i) = tors6(1,i) * tlambda
//                   end if
//                end do
//             end if
//          end do
//       end if
//       return
//       end
