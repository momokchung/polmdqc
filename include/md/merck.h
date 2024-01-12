// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  merck  --  MMFF-specific force field parameters  //
//                                                   //
///////////////////////////////////////////////////////

// nligne     number of atom pairs having MMFF Bond Type 1
// bt_1       atom pairs having MMFF Bond Type 1
// eqclass    table of atom class equivalencies used to find
//            default parameters if explicit values are missing
//            (see J. Comput. Chem., 17, 490-519, '95, Table IV)
// crd        number of attached neighbors    |
// val        valency value                   |  see T. A. Halgren,
// pilp       if 0, no lone pair              |  J. Comput. Chem.,
//            if 1, one or more lone pair(s)  |  17, 616-645 (1995)
// mltb       multibond indicator             |
// arom       aromaticity indicator           |
// lin        linearity indicator             |
// sbmb       single- vs multiple-bond flag   |
// mmffarom   aromatic rings parameters
// mmffaromc  cationic aromatic rings parameters
// mmffaroma  anionic aromatic rings parameters

MDQC_EXTERN int nligne ;
MDQC_EXTERN int bt_1[2][500];
MDQC_EXTERN int eqclass[5][500];
MDQC_EXTERN int crd[100];
MDQC_EXTERN int val[100];
MDQC_EXTERN int pilp[100];
MDQC_EXTERN int mltb[100];
MDQC_EXTERN int arom[100];
MDQC_EXTERN int lin[100];
MDQC_EXTERN int sbmb[100];
MDQC_EXTERN int mmffarom[6][maxtyp];
MDQC_EXTERN int mmffaromc[6][maxtyp];
MDQC_EXTERN int mmffaroma[6][maxtyp];

// rad0      covalent atomic radius for empirical bond rules
// paulel    Pauling electronegativities for empirical bond rules
// r0ref     reference bond length for empirical bond rules
// kbref     reference force constant for empirical bond rules
// mmff_kb   bond force constant for pairs of atom classes
// mmff_kb1  bond force constant for class pairs with Bond Type 1
// mmff_b0   bond length value for pairs of atom classes
// mmff_b1   bond length value for class pairs with Bond Type 1

MDQC_EXTERN real rad0[100];
MDQC_EXTERN real paulel[100];
MDQC_EXTERN real r0ref[100][100];
MDQC_EXTERN real kbref[100][100];
MDQC_EXTERN real mmff_kb[100][100];
MDQC_EXTERN real mmff_kb1[100][100];
MDQC_EXTERN real mmff_b0[100][100];
MDQC_EXTERN real mmff_b1[100][100];

// mmff_ka     angle force constant for triples of atom classes
// mmff_ka1    angle force constant with one bond of Type 1
// mmff_ka2    angle force constant with both bonds of Type 1
// mmff_ka3    angle force constant for 3-membered ring
// mmff_ka4    angle force constant for 4-membered ring
// mmff_ka5    angle force constant for 3-ring and one Bond Type 1
// mmff_ka6    angle force constant for 3-ring and both Bond Type 1
// mmff_ka7    angle force constant for 4-ring and one Bond Type 1
// mmff_ka8    angle force constant for 4-ring and both Bond Type 1
// mmff_ang0   ideal bond angle for triples of atom classes
// mmff_ang1   ideal bond angle with one bond of Type 1
// mmff_ang2   ideal bond angle with both bonds of Type 1
// mmff_ang3   ideal bond angle for 3-membered ring
// mmff_ang4   ideal bond angle for 4-membered ring
// mmff_ang5   ideal bond angle for 3-ring and one Bond Type 1
// mmff_ang6   ideal bond angle for 3-ring and both Bond Type 1
// mmff_ang7   ideal bond angle for 4-ring and one Bond Type 1
// mmff_ang8   ideal bond angle for 4-ring and both Bond Type 1

MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka1;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka2;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka3;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka4;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka5;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka6;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka7;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ka8;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang0;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang1;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang2;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang3;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang4;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang5;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang6;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang7;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> mmff_ang8;

// Stretch-Bend Type 0
// stbn_abc     stretch-bend parameters for A-B-C atom classes
// stbn_cba     stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 1  (A-B is Bond Type 1)
// stbn_abc1    stretch-bend parameters for A-B-C atom classes
// stbn_cba1    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 2  (B-C is Bond Type 1) 
// stbn_abc2    stretch-bend parameters for A-B-C atom classes
// stbn_cba2    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type = 3  (A-B and B-C are Bond Type 1) 
// stbn_abc3    stretch-bend parameters for A-B-C atom classes
// stbn_cba3    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 4  (both Bond Types 0, 4-membered ring)
// stbn_abc4    stretch-bend parameters for A-B-C atom classes
// stbn_cba4    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 5  (both Bond Types 0, 3-membered ring)
// stbn_abc5    stretch-bend parameters for A-B-C atom classes
// stbn_cba5    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 6  (A-B is Bond Type 1, 3-membered ring)
// stbn_abc6    stretch-bend parameters for A-B-C atom classes
// stbn_cba6    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 7  (B-C is Bond Type 1, 3-membered ring)
// stbn_abc7    stretch-bend parameters for A-B-C atom classes
// stbn_cba7    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 8  (both Bond Types 1, 3-membered ring)
// stbn_abc8    stretch-bend parameters for A-B-C atom classes
// stbn_cba8    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 9  (A-B is Bond Type 1, 4-membered ring)
// stbn_abc9    stretch-bend parameters for A-B-C atom classes
// stbn_cba9    stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 10  (B-C is Bond Type 1, 4-membered ring)
// stbn_abc10   stretch-bend parameters for A-B-C atom classes
// stbn_cba10   stretch-bend parameters for C-B-A atom classes
// Stretch-Bend Type 11  (both Bond Types 1, 4-membered ring)
// stbn_abc11   stretch-bend parameters for A-B-C atom classes
// stbn_cba11   stretch-bend parameters for C-B-A atom classes
// defstbn_abc  default stretch-bend parameters for A-B-C classes
// defstbn_cba  default stretch-bend parameters for C-B-A classes

MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc1;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba1;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc2;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba2;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc3;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba3;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc4;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba4;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc5;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba5;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc6;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba6;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc7;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba7;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc8;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba8;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc9;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba9;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc10;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba10;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_abc11;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> stbn_cba11;
MDQC_EXTERN real defstbn_abc[5][5][5];
MDQC_EXTERN real defstbn_cba[5][5][5];

// t1_1     torsional parameters for 1-fold, MMFF Torsion Type 1
// t1_2     torsional parameters for 1-fold, MMFF Torsion Type 2
// t2_1     torsional parameters for 2-fold, MMFF Torsion Type 1
// t2_2     torsional parameters for 2-fold, MMFF Torsion Type 2
// t3_1     torsional parameters for 3-fold, MMFF Torsion Type 1
// t3_2     torsional parameters for 3-fold, MMFF Torsion Type 2
// kt_1     string of classes for torsions, MMFF Torsion Type 1
// kt_2     string of classes for torsions, MMFF Torsion Type 2

MDQC_EXTERN real t1_1[2001][2];
MDQC_EXTERN real t2_1[2001][2];
MDQC_EXTERN real t3_1[2001][2];
MDQC_EXTERN real t1_2[2001][2];
MDQC_EXTERN real t2_2[2001][2];
MDQC_EXTERN real t3_2[2001][2];
MDQC_EXTERN std::string kt_1[2001];
MDQC_EXTERN std::string kt_2[2001];

// g        scale factors for calculation of MMFF eps
// alph     atomic polarizabilities for calculation of MMFF eps
// nn       effective number of valence electrons for MMFF eps
// da       donor/acceptor atom classes

MDQC_EXTERN real g[maxclass];
MDQC_EXTERN real alph[maxclass];
MDQC_EXTERN real nn[maxclass];
MDQC_EXTERN char da[maxclass];

// bci      bond charge increments for building atom charges
// bci_1    bond charge increments for MMFF Bond Type 1   
// pbci     partial BCI for building missing BCI's
// fcadj    formal charge adjustment factor

MDQC_EXTERN real bci[100][100];
MDQC_EXTERN real bci_1[100][100];
MDQC_EXTERN real pbci[maxclass];
MDQC_EXTERN real fcadj[maxclass];
}
