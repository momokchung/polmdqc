// Author: Moses KJ Chung
// Year:   2024

#include "basis.h"
#include "gettext.h"
#include "kgbs.h"
#include "ptable.h"
#include "readbasis.h"
#include "upcase.h"
#include <sstream>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  readbasis  --  input of basis set parameters  //
//                                                //
////////////////////////////////////////////////////

// "readbasis" processes the basis set parameter file

void readbasis()
{
    int next;
    int ibss;
    int atmind,nprim,nbasis,am;
    realq scale,exp,coeff1,coeff2;
    BasisType bsstyp;
    std::string keyword,record,string;
    std::string elename,ams;
    std::istringstream iss;

    // set default bsstyp
    bsstyp = BasisType::Spherical;

    // initialize number of basis functions
    for (int i = 0; i < maxele; i++) {
        ngbs[i] = 0;
    }

    // process each line of the basis set file
    ibss = 0;
    while (ibss < nbss) {
        record = bssline[ibss];
        next = 0;
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);

        if (keyword == "SPHERICAL") {
            bsstyp = BasisType::Spherical;
        }
        else if (keyword == "CARTESIAN") {
            bsstyp = BasisType::Cartesian;
        }
        else if (keyword == "****") {
            ibss++;
            record = bssline[ibss];
            iss.clear();
            iss.str(record);

            // get element name
            elename = "";
            iss >> elename;
            upcase(elename);
            if (elename == "") continue;
            atmind = symtoatmn[elename] - 1;

            // read basis set
            keyword = "";
            nbasis = 0;
            while (keyword != "****" and ibss < nbss) {
                ibss++;
                record = bssline[ibss];
                iss.clear();
                iss.str(record);

                // get angular momentum, # prim, scale
                ams = "";
                nprim = 0;
                scale = 0.;
                iss >> ams >> nprim >> scale;
                upcase(ams);
                if (ams == "") continue;

                // read basis function
                if (ams == "SP") {
                    amgbs[atmind][nbasis] = 0;
                    nprimgbs[atmind][nbasis] = nprim;
                    scalegbs[atmind][nbasis] = scale;
                    amgbs[atmind][nbasis+1] = 1;
                    nprimgbs[atmind][nbasis+1] = nprim;
                    scalegbs[atmind][nbasis+1] = scale;

                    // read primitive function
                    for (int i = 0; i < nprim; i++) {
                        ibss++;
                        record = bssline[ibss];
                        iss.clear();
                        iss.str(record);
                        iss >> exp >> coeff1 >> coeff2;
                        expgbs[atmind][nbasis][i] = exp;
                        expgbs[atmind][nbasis+1][i] = exp;
                        coeffgbs[atmind][nbasis][i] = coeff1;
                        coeffgbs[atmind][nbasis+1][i] = coeff2;
                    }
                    nbasis += 2;
                }
                else {
                    if (ams == "S") am = 0;
                    else if (ams == "P") am = 1;
                    else if (ams == "D") am = 2;
                    else if (ams == "F") am = 3;
                    else if (ams == "G") am = 4;
                    else if (ams == "H") am = 5;
                    else if (ams == "I") am = 6;
                    else if (ams == "K" or ams == "L=7") am = 7;
                    amgbs[atmind][nbasis] = am;
                    nprimgbs[atmind][nbasis] = nprim;
                    scalegbs[atmind][nbasis] = scale;

                    // read primitive function
                    for (int i = 0; i < nprim; i++) {
                        ibss++;
                        record = bssline[ibss];
                        iss.clear();
                        iss.str(record);
                        iss >> exp >> coeff1;
                        expgbs[atmind][nbasis][i] = exp;
                        coeffgbs[atmind][nbasis][i] = coeff1;
                    }
                    nbasis++;
                }

                // peak if next line is "****"
                record = bssline[ibss+1];
                iss.clear();
                iss.str(record);
                iss >> keyword;
            }
            ngbs[atmind] = nbasis;
        }
        ibss++;
    }

    // define shared parameters for each element
    for (int i = 0; i < maxele; i++) {
        if (ngbs[i] > 0) {
            namegbs[i] = bssname;
            typgbs[i] = bsstyp;
        }
    }

    // // uncomment below to check readgbs
    // for (int i = 0; i < maxele; i++) {
    //     if (ngbs[i] > 0) {
    //         printf("\n Element : %s", atmntosym[i+1].c_str());
    //         printf("\n # Basis Functions : %d", ngbs[i]);
    //         for (int j = 0; j < ngbs[i]; j++) {
    //             printf("\n    AngMom : %d", amgbs[i][j]);
    //             printf("\n     Scale : %10.6f", scalegbs[i][j]);
    //             printf("\n    # Prim : %d", nprimgbs[i][j]);
    //             printf("\n         Exp : ");
    //             for (int k = 0; k < nprimgbs[i][j]; k++) {
    //                 printf("%15.7f", expgbs[i][j][k]);
    //             }
    //             printf("\n       Coeff : ");
    //             for (int k = 0; k < nprimgbs[i][j]; k++) {
    //                 printf("%15.7f", coeffgbs[i][j][k]);
    //             }
    //             printf("\n");
    //         }
    //         printf("\n****\n");
    //     }
    // }
}
}
