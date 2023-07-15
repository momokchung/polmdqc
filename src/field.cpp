/////////////////////////////////////////////////////////
//                                                     //
//  field.cpp  --  get the potential energy functions  //
//                                                     //
/////////////////////////////////////////////////////////

// "field" sets the force field potential energy functions from
// a parameter file and modifications specified in a keyfile


#include "field.h"
#include "fields.h"
#include "getprm.h"
#include "inform.h"
#include "keys.h"
#include "potent.h"
#include "sizes.h"
#include <string>

void field()
{
    std::string keyword;
    std::string record;
    std::string string;
    int next;
    bool header;

    // set the default values for the active potentials
    use_bond = true;
    use_angle = true;
    use_strbnd = true;
    use_urey = true;
    use_angang = true;
    use_opbend = true;
    use_opdist = true;
    use_improp = true;
    use_imptor = true;
    use_tors = true;
    use_pitors = true;
    use_strtor = true;
    use_angtor = true;
    use_tortor = true;
    use_vdw = true;
    use_repel = true;
    use_disp = true;
    use_charge = true;
    use_chgdpl = true;
    use_dipole = true;
    use_mpole = true;
    use_polar = true;
    use_chgtrn = true;
    use_chgflx = true;
    use_rxnfld = false;
    use_solv = true;
    use_metal = false;
    use_geom = true;
    use_extra = true;

    // read the potential energy force field parameter file
    getprm();
}



// c
// c     check keywords for biopolymer atom type definitions
// c
//       header = .true.
//       do i = 1, nkey
//          next = 1
//          record = keyline(i)
//          call gettext (record,keyword,next)
//          call upcase (keyword)
//          if (keyword(1:8) .eq. 'BIOTYPE ') then
//             ia = 0
//             ib = 0
//             string = record(next:240)
//             read (string,*,err=10,end=10)  ia
//             call getword (record,string,next)
//             call getstring (record,string,next)
//             string = record(next:240)
//             read (string,*,err=10,end=10)  ib
//    10       continue
//             if (ia.ge.0 .and. ia.le.maxbio) then
//                if (header .and. .not.silent) then
//                   header = .false.
//                   write (iout,20)
//    20             format (/,' Additional Biopolymer Type Definitions',
//      &                    //,5x,'BioType',10x,'Atom Type',/)
//                end if
//                biotyp(ia) = ib
//                if (.not. silent) then
//                   write (iout,30)  ia,ib
//    30             format (1x,i8,8x,i8)
//                end if

//             else if (ia .gt. maxbio) then
//                write (iout,40)
//    40          format (/,' FIELDS  --  Too many Biopolymer Types;',
//      &                    ' Increase MAXBIO')
//                call fatal
//             end if
//          end if
//       end do
// c
// c     check keywords for potential function control parameters
// c
//       do i = 1, nkey
//          record = keyline(i)
//          call prmkey (record)
//       end do
//       return
//       end
