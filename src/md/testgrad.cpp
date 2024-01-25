// Author: Moses KJ Chung
// Year:   2024

#include "atoms.h"
#include "calcMode.h"
#include "deriv.h"
#include "energi.h"
#include "energy.h"
#include "files.h"
#include "final.h"
#include "getcart.h"
#include "gettext.h"
#include "inform.h"
#include "initial.h"
#include "libfunc.h"
#include "mechanic.h"
#include "nextarg.h"
#include "output.h"
#include "readcart.h"
#include "upcase.h"
#include "usage.h"
#include <iostream>
#include <sstream>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////
//                                 //
//  testgrad  --  derivative test  //
//                                 //
/////////////////////////////////////

// "testgrad" computes and compares the analytical and numerical
// gradient vectors of the potential energy function with respect
// to Cartesian coordinates

void resizeNumDer();

void testgrad(int argc, char** argv)
{
    constexpr CalcMode EnergyMode = CalcMode::Energy;
    constexpr CalcMode GradientMode = CalcMode::Gradient;

    // integer i,j,ixyz
    int next,frame;
    int nold;
    int numSpaces;
    int width;
    int precision;
    // integer freeunit
    real f,f0,eps,eps0,old;
    real eb0,ea0,eba0,eub0,eaa0,eopb0;
    real eopd0,eid0,eit0,et0,ept0,ebt0;
    real eat0,ett0,ev0,er0,edsp0,ec0;
    real ecd0,ed0,em0,ep0,ect0,erxf0;
    real es0,elf0,eg0,ex0;
    real totnorm,ntotnorm,rms,nrms;
    std::vector<real> denorm;
    std::vector<real> ndenorm;
    std::vector<real> told;
    bool exist,query;
    bool doanalyt,donumer,dofull;
    bool first;
    // character*1 axis(3)
    std::string answer;
    std::string record,string;
    std::istringstream iss;
    const char axis[] = {'X', 'Y', 'Z'};

    // set up the structure and mechanics calculation
    initial(argc, argv);
    getcart(ffile);
    mechanic();

    // decide whether to do an analytical gradient calculation
    doanalyt = true;
    nextarg(answer,exist);
    if (!exist) {
        printf("\n Compute the Analytical Gradient Vector [Y] :  ");
        std::getline(std::cin, record);
        next = 0;
        gettext(record,answer,next);
    }
    upcase(answer);
    if (answer == "N") doanalyt = false;

    // decide whether to do a numerical gradient calculation
    donumer = true;
    nextarg(answer,exist);
    if (!exist) {
        printf("\n Compute the Numerical Gradient Vector [Y] :   ");
        std::getline(std::cin, record);
        next = 0;
        gettext(record,answer,next);
    }
    upcase(answer);
    if (answer == "N") donumer = false;

    // get the stepsize for numerical gradient calculation
    if (donumer) {
        eps = -1.;
        eps0 = 0.00001;
        query = true;
        nextarg(string,exist);
        if (exist) {
            iss.clear();
            iss.str(string);
            iss >> eps;
            query = false;
        }
        if (query) {
            printf("\n Enter Finite Difference Stepsize [%8.1e Ang] :  ", eps0);
            std::getline(std::cin, record);
            iss.clear();
            iss.str(record);
            iss >> eps;
        }
        if (eps < 0.) eps = eps0;
    }

    // decide whether to output results by gradient component
    dofull = true;
    if (n > 100) {
        dofull = false;
        nextarg(answer,exist);
        if (!exist) {
            printf("\n Output Breakdown by Gradient Component [N] :  ");
            std::getline(std::cin, record);
            next = 0;
            gettext(record,answer,next);
        }
        upcase(answer);
        if (answer == "Y") dofull = true;
    }

    // perform dynamic allocation of some local arrays
    denorm.resize(n);
    if (donumer) {
        ndenorm.resize(n);
        resizeNumDer();
    }

    // perform analysis for each successive coordinate structure
    frame = 0;
    while (!informAbort) {
        frame++;
        if (frame > 1) {
            if (!test) printf("\n Analysis for Archive Structure :        %8d\n", frame);
            if (nold != n) {
                mechanic();
                denorm.resize(n);
                if (donumer) {
                    ndenorm.resize(n);
                    resizeNumDer();
                }
            }
            else {
                for (int i = 0; i < n; i++) {
                    if (type[i] != told[i]) {
                        mechanic();
                        break;                        
                    }
                }
            }
        }

        // compute the analytical gradient components
        if (doanalyt) {
            energy<GradientMode>();
        }

        // print the total potential energy of the system
        if (doanalyt and !test) {
            numSpaces = 8;
            if (digits >= 8) {
                width = 20;
                precision = 8;
            }
            else if (digits >= 6) {
                width = 18;
                precision = 6;
            }
            else {
                width = 16;
                precision = 4;
            }
            printf("\n Total Potential Energy :%*s%*.*f Kcal/mole\n", numSpaces, "", width, precision, esum);


            // print the energy breakdown over individual components
            printf("\n Potential Energy Breakdown by Individual Components :\n");

            if (digits >= 8) {
                printf("\n  Energy       EB              EA              EBA             EUB");
                printf("\n  Terms        EAA             EOPB            EOPD            EID");
                printf("\n               EIT             ET              EPT             EBT");
                printf("\n               EAT             ETT             EV              ER");
                printf("\n               EDSP            EC              ECD             ED");
                printf("\n               EM              EP              ECT             ERXF");
                printf("\n               ES              ELF             EG              EX\n");
                printf("\n      %16.8f%16.8f%16.8f%16.8f", eb, ea, eba, eub);
                printf("\n      %16.8f%16.8f%16.8f%16.8f", eaa, eopb, eopd, eid);
                printf("\n      %16.8f%16.8f%16.8f%16.8f", eit, et, ept, ebt);
                printf("\n      %16.8f%16.8f%16.8f%16.8f", eat, ett, ev, er);
                printf("\n      %16.8f%16.8f%16.8f%16.8f", edsp, ec, ecd, ed);
                printf("\n      %16.8f%16.8f%16.8f%16.8f", em, ep, ect, erxf);
                printf("\n      %16.8f%16.8f%16.8f%16.8f\n", es, elf, eg, ex);
            }
            else if (digits >= 6) {
                printf("\n  Energy      EB            EA            EBA           EUB           EAA");
                printf("\n  Terms       EOPB          EOPD          EID           EIT           ET");
                printf("\n              EPT           EBT           EAT           ETT           EV");
                printf("\n              ER            EDSP          EC            ECD           ED");
                printf("\n              EM            EP            ECT           ERXF          ES");
                printf("\n              ELF           EG            EX\n");
                printf("\n      %14.6f%14.6f%14.6f%14.6f%14.6f", eb, ea, eba, eub, eaa);
                printf("\n      %14.6f%14.6f%14.6f%14.6f%14.6f", eopb, eopd, eid, eit, et);
                printf("\n      %14.6f%14.6f%14.6f%14.6f%14.6f", ept, ebt, eat, ett, ev);
                printf("\n      %14.6f%14.6f%14.6f%14.6f%14.6f", er, edsp, ec, ecd, ed);
                printf("\n      %14.6f%14.6f%14.6f%14.6f%14.6f", em, ep, ect, erxf, es);
                printf("\n      %14.6f%14.6f%14.6f\n", elf, eg, ex);
            }
            else {
                printf("\n  Energy      EB          EA          EBA         EUB         EAA         EOPB");
                printf("\n  Terms       EOPD        EID         EIT         ET          EPT         EBT");
                printf("\n              EAT         ETT         EV          ER          EDSP        EC");
                printf("\n              ECD         ED          EM          EP          ECT         ERXF");
                printf("\n              ES          ELF         EG          EX\n");
                printf("\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f", eb, ea, eba, eub, eaa, eopb);
                printf("\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f", eopd, eid, eit, et, ept, ebt);
                printf("\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f", eat, ett, ev, er, edsp, ec);
                printf("\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f", ecd, ed, em, ep, ect, erxf);
                printf("\n      %12.4f%12.4f%12.4f%12.4f\n", es, elf, eg, ex);
            }
        }

        // print a header for the gradients of individual potentials
        if (dofull and !test) {
            printf("\n Cartesian Gradient Breakdown by Individual Components :\n");
            if (digits >= 8) {
                printf("\n  Atom         d EB            d EA            d EBA           d EUB");
                printf("\n  Axis         d EAA           d EOPB          d EOPD          d EID");
                printf("\n  Type         d EIT           d ET            d EPT           d EBT");
                printf("\n               d EAT           d ETT           d EV            d ER");
                printf("\n               d EDSP          d EC            d ECD           d ED");
                printf("\n               d EM            d EP            d ECT           d ERXF");
                printf("\n               d ES            d ELF           d EG            d EX\n");
            }
            else if (digits >= 6) {
                    printf("\n  Atom        d EB          d EA          d EBA         d EUB         d EAA");
                    printf("\n  Axis        d EOPB        d EOPD        d EID         d EIT         d ET");
                    printf("\n  Type        d EPT         d EBT         d EAT         d ETT         d EV");
                    printf("\n              d ER          d EDSP        d EC          d ECD         d ED");
                    printf("\n              d EM          d EP          d ECT         d ERXF        d ES");
                    printf("\n              d ELF         d EG          d EX\n");
            }
            else {
                printf("\n  Atom      d EB        d EA        d EBA       d EUB       d EAA       d EOPB");
                printf("\n  Axis      d EOPD      d EID       d EIT       d ET        d EPT       d EBT");
                printf("\n  Type      d EAT       d ETT       d EV        d ER        d EDSP      d EC");
                printf("\n            d ECD       d ED        d EM        d EP        d ECT       d ERXF");
                printf("\n            d ES        d ELF       d EG        d EX\n");
            }
        }

        // get the Cartesian component two-sided numerical gradients
        for (int i = 0; i < n; i++) {
            if (donumer and use[i+1]) {
                for (int j = 0; j < 3; j++) {
                    if (j == 0) {
                        old = x[i];
                        x[i] -= (real)0.5 * eps;
                    }
                    else if (j == 1) {
                        old = y[i];
                        y[i] -= (real)0.5 * eps;
                    }
                    else if (j == 2) {
                        old = z[i];
                        z[i] -= (real)0.5 * eps;
                    }
                    energy<EnergyMode>();
                    f0 = esum;
                    eb0 = eb;
                    ea0 = ea;
                    eba0 = eba;
                    eub0 = eub;
                    eaa0 = eaa;
                    eopb0 = eopb;
                    eopd0 = eopd;
                    eid0 = eid;
                    eit0 = eit;
                    et0 = et;
                    ept0 = ept;
                    ebt0 = ebt;
                    eat0 = eat;
                    ett0 = ett;
                    ev0 = ev;
                    er0 = er;
                    edsp0 = edsp;
                    ec0 = ec;
                    ecd0 = ecd;
                    ed0 = ed;
                    em0 = em;
                    ep0 = ep;
                    ect0 = ect;
                    erxf0 = erxf;
                    es0 = es;
                    elf0 = elf;
                    eg0 = eg;
                    ex0 = ex;
                    if (j == 0) {
                        x[i] += eps;
                    }
                    else if (j == 1) {
                        y[i] += eps;
                    }
                    else if (j == 2) {
                        z[i] += eps;
                    }
                    energy<EnergyMode>();
                    f = esum;
                    if (j == 0) {
                        x[i] = old;
                    }
                    else if (j == 1) {
                        y[i] = old;
                    }
                    else if (j == 2) {
                        z[i] = old;
                    }
                    ndesum[3*i+j] = (f - f0) / eps;
                    ndeb[3*i+j] = (eb - eb0) / eps;
                    ndea[3*i+j] = (ea - ea0) / eps;
                    ndeba[3*i+j] = (eba - eba0) / eps;
                    ndeub[3*i+j] = (eub - eub0) / eps;
                    ndeaa[3*i+j] = (eaa - eaa0) / eps;
                    ndeopb[3*i+j] = (eopb - eopb0) / eps;
                    ndeopd[3*i+j] = (eopd - eopd0) / eps;
                    ndeid[3*i+j] = (eid - eid0) / eps;
                    ndeit[3*i+j] = (eit - eit0) / eps;
                    ndet[3*i+j] = (et - et0) / eps;
                    ndept[3*i+j] = (ept - ept0) / eps;
                    ndebt[3*i+j] = (ebt - ebt0) / eps;
                    ndeat[3*i+j] = (eat - eat0) / eps;
                    ndett[3*i+j] = (ett - ett0) / eps;
                    ndev[3*i+j] = (ev - ev0) / eps;
                    nder[3*i+j] = (er - er0) / eps;
                    ndedsp[3*i+j] = (edsp - edsp0) / eps;
                    ndec[3*i+j] = (ec - ec0) / eps;
                    ndecd[3*i+j] = (ecd - ecd0) / eps;
                    nded[3*i+j] = (ed - ed0) / eps;
                    ndem[3*i+j] = (em - em0) / eps;
                    ndep[3*i+j] = (ep - ep0) / eps;
                    ndect[3*i+j] = (ect - ect0) / eps;
                    nderxf[3*i+j] = (erxf - erxf0) / eps;
                    ndes[3*i+j] = (es - es0) / eps;
                    ndelf[3*i+j] = (elf - elf0) / eps;
                    ndeg[3*i+j] = (eg - eg0) / eps;
                    ndex[3*i+j] = (ex - ex0) / eps;
                }
            }

            // print analytical gradients of each energy term for each atom
            if (dofull and use[i+1] and !test) {
                for (int j = 0; j < 3; j++) {
                    if (doanalyt) {
                        if (digits >= 8) {
                            printf("\n%6d%16.8f%16.8f%16.8f%16.8f"
                            "\n     %c%16.8f%16.8f%16.8f%16.8f"
                                "\n Anlyt%16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f\n",
                                i+1,deb[3*i+j],dea[3*i+j],deba[3*i+j],deub[3*i+j],
                                axis[j],deaa[3*i+j],deopb[3*i+j],deopd[3*i+j],deid[3*i+j],
                                deit[3*i+j],det[3*i+j],dept[3*i+j],debt[3*i+j],
                                deat[3*i+j],dett[3*i+j],dev[3*i+j],der[3*i+j],
                                dedsp[3*i+j],dec[3*i+j],decd[3*i+j],ded[3*i+j],
                                dem[3*i+j],dep[3*i+j],dect[3*i+j],derxf[3*i+j],
                                des[3*i+j],delf[3*i+j],deg[3*i+j],dex[3*i+j]);
                        }
                        else if (digits >= 6) {
                            printf("\n%6d%14.6f%14.6f%14.6f%14.6f%14.6f"
                            "\n     %c%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n Anlyt%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f\n",
                                i+1,deb[3*i+j],dea[3*i+j],deba[3*i+j],deub[3*i+j],deaa[3*i+j],
                                axis[j],deopb[3*i+j],deopd[3*i+j],deid[3*i+j],deit[3*i+j],det[3*i+j],
                                dept[3*i+j],debt[3*i+j],deat[3*i+j],dett[3*i+j],dev[3*i+j],
                                der[3*i+j],dedsp[3*i+j],dec[3*i+j],decd[3*i+j],ded[3*i+j],
                                dem[3*i+j],dep[3*i+j],dect[3*i+j],derxf[3*i+j],des[3*i+j],
                                delf[3*i+j],deg[3*i+j],dex[3*i+j]);
                        }
                        else {
                            printf("\n%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                            "\n     %c%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n Anlyt%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f\n",
                                i+1,deb[3*i+j],dea[3*i+j],deba[3*i+j],deub[3*i+j],deaa[3*i+j],deopb[3*i+j],
                                axis[j],deopd[3*i+j],deid[3*i+j],deit[3*i+j],det[3*i+j],dept[3*i+j],debt[3*i+j],
                                deat[3*i+j],dett[3*i+j],dev[3*i+j],der[3*i+j],dedsp[3*i+j],dec[3*i+j],
                                decd[3*i+j],ded[3*i+j],dem[3*i+j],dep[3*i+j],dect[3*i+j],derxf[3*i+j],
                                des[3*i+j],delf[3*i+j],deg[3*i+j],dex[3*i+j]);
                        }
                    }

                    // print numerical gradients of each energy term for each atom
                    if (donumer) {
                        if (digits >= 8) {
                            printf("\n%6d%16.8f%16.8f%16.8f%16.8f"
                            "\n     %c%16.8f%16.8f%16.8f%16.8f"
                                "\n Numer%16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f"
                                "\n      %16.8f%16.8f%16.8f%16.8f\n",
                                i+1,ndeb[3*i+j],ndea[3*i+j],ndeba[3*i+j],ndeub[3*i+j],
                                axis[j],ndeaa[3*i+j],ndeopb[3*i+j],ndeopd[3*i+j],ndeid[3*i+j],
                                ndeit[3*i+j],ndet[3*i+j],ndept[3*i+j],ndebt[3*i+j],
                                ndeat[3*i+j],ndett[3*i+j],ndev[3*i+j],nder[3*i+j],
                                ndedsp[3*i+j],ndec[3*i+j],ndecd[3*i+j],nded[3*i+j],
                                ndem[3*i+j],ndep[3*i+j],ndect[3*i+j],nderxf[3*i+j],
                                ndes[3*i+j],ndelf[3*i+j],ndeg[3*i+j],ndex[3*i+j]);
                        }
                        else if (digits >= 6) {
                            printf("\n%6d%14.6f%14.6f%14.6f%14.6f%14.6f"
                            "\n     %c%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n Numer%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f\n",
                                i+1,ndeb[3*i+j],ndea[3*i+j],ndeba[3*i+j],ndeub[3*i+j],ndeaa[3*i+j],
                                axis[j],ndeopb[3*i+j],ndeopd[3*i+j],ndeid[3*i+j],ndeit[3*i+j],ndet[3*i+j],
                                ndept[3*i+j],ndebt[3*i+j],ndeat[3*i+j],ndett[3*i+j],ndev[3*i+j],
                                nder[3*i+j],ndedsp[3*i+j],ndec[3*i+j],ndecd[3*i+j],nded[3*i+j],
                                ndem[3*i+j],ndep[3*i+j],ndect[3*i+j],nderxf[3*i+j],ndes[3*i+j],
                                ndelf[3*i+j],ndeg[3*i+j],ndex[3*i+j]);
                        }
                        else {
                            printf("\n%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                            "\n     %c%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n Numer%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f\n",
                                i+1,ndeb[3*i+j],ndea[3*i+j],ndeba[3*i+j],ndeub[3*i+j],ndeaa[3*i+j],ndeopb[3*i+j],
                                axis[j],ndeopd[3*i+j],ndeid[3*i+j],ndeit[3*i+j],ndet[3*i+j],ndept[3*i+j],ndebt[3*i+j],
                                ndeat[3*i+j],ndett[3*i+j],ndev[3*i+j],nder[3*i+j],ndedsp[3*i+j],ndec[3*i+j],
                                ndecd[3*i+j],nded[3*i+j],ndem[3*i+j],ndep[3*i+j],ndect[3*i+j],nderxf[3*i+j],
                                ndes[3*i+j],ndelf[3*i+j],ndeg[3*i+j],ndex[3*i+j]);
                        }
                    }
                }
            }
        }

        // print the total gradient components for each atom
        if ((doanalyt or donumer) and !test) {
            printf("\n Cartesian Gradient Breakdown over Individual Atoms :\n");
            int s1,s2,s3,s4,s5;
            if (digits >= 8) {
                s1 = 4; s2 = 10; s3 = 11; s4 = 11; s5 = 11;
            }
            else if (digits >= 6) {
                s1 = 6; s2 = 11; s3 = 9; s4 = 9; s5 = 11;
            }
            else {
                s1 = 6; s2 = 14; s3 = 7; s4 = 7; s5 = 10;
            }
            printf("\n  Type%*sAtom%*sdE/dX%*sdE/dY%*sdE/dZ%*sNorm\n\n", s1,"",s2,"",s3,"",s4,"",s5,"");
        }
        if (!test) {
            totnorm = 0.;
            ntotnorm = 0.;
            for (int i = 0; i < n; i++) {
                if (doanalyt and use[i+1]) {
                    denorm[i] = REAL_POW(desum[3*i+0],2) + REAL_POW(desum[3*i+1],2) + REAL_POW(desum[3*i+2],2);
                    totnorm = totnorm + denorm[i];
                    denorm[i] = REAL_SQRT(denorm[i]);
                    if (digits >= 8) {
                        printf(" Anlyt%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,desum[3*i+0],desum[3*i+1],desum[3*i+2],denorm[i]);
                    }
                    else if (digits >= 6) {
                        printf(" Anlyt  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,desum[3*i+0],desum[3*i+1],desum[3*i+2],denorm[i]);
                    }
                    else {
                        printf(" Anlyt  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,desum[3*i+0],desum[3*i+1],desum[3*i+2],denorm[i]);
                    }
                }
                if (donumer and use[i+1]) {
                    ndenorm[i] = REAL_POW(ndesum[3*i+0],2) + REAL_POW(ndesum[3*i+1],2) + REAL_POW(ndesum[3*i+2],2);
                    ntotnorm = ntotnorm + ndenorm[i];
                    ndenorm[i] = REAL_SQRT(ndenorm[i]);
                    if (digits >= 8) {
                        printf(" Numer%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,ndesum[3*i+0],ndesum[3*i+1],ndesum[3*i+2],ndenorm[i]);
                    }
                    else if (digits >= 6) {
                        printf(" Numer  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,ndesum[3*i+0],ndesum[3*i+1],ndesum[3*i+2],ndenorm[i]);
                    }
                    else {
                        printf(" Numer  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,ndesum[3*i+0],ndesum[3*i+1],ndesum[3*i+2],ndenorm[i]);
                    }
                }
            }

            // print the total norm for the analytical gradient
            if (doanalyt or donumer) {
                printf("\n Total Gradient Norm and RMS Gradient per Atom :\n\n");
            }
            if (doanalyt) {
                totnorm = REAL_SQRT(totnorm);
                if (digits >= 8) {
                    printf(" Anlyt      Total Gradient Norm Value      %20.8f\n", totnorm);
                }
                else if (digits >= 6) {
                    printf(" Anlyt      Total Gradient Norm Value      %18.6f\n", totnorm);
                }
                else {
                    printf(" Anlyt      Total Gradient Norm Value      %16.4f\n", totnorm);
                }
            }

            // print the total norm for the numerical gradient
            if (donumer) {
                ntotnorm = REAL_SQRT(ntotnorm);
                if (digits >= 8) {
                    printf(" Numer      Total Gradient Norm Value      %20.8f\n", ntotnorm);
                }
                else if (digits >= 6) {
                    printf(" Numer      Total Gradient Norm Value      %18.6f\n", ntotnorm);
                }
                else {
                    printf(" Numer      Total Gradient Norm Value      %16.4f\n", ntotnorm);
                }
            }

            printf("\n");

            // print the rms per atom norm for the analytical gradient
            if (doanalyt) {
                rms = totnorm / REAL_SQRT(static_cast<real>(nuse));
                if (digits >= 8) {
                    printf(" Anlyt      RMS Gradient over All Atoms    %20.8f\n", rms);
                }
                else if (digits >= 6) {
                    printf(" Anlyt      RMS Gradient over All Atoms    %18.6f\n", rms);
                }
                else {
                    printf(" Anlyt      RMS Gradient over All Atoms    %16.4f\n", rms);
                }
            }

            // print the rms per atom norm for the numerical gradient
            if (donumer) {
                nrms = ntotnorm / REAL_SQRT(static_cast<real>(nuse));
                if (digits >= 8) {
                    printf(" Numer      RMS Gradient over All Atoms    %20.8f\n", nrms);
                }
                else if (digits >= 6) {
                    printf(" Numer      RMS Gradient over All Atoms    %18.6f\n", nrms);
                }
                else {
                    printf(" Numer      RMS Gradient over All Atoms    %16.4f\n", nrms);
                }
            }
        }

        // attempt to read next structure from the coordinate file
        if (told.size() != n) {
            told.resize(n);
        }
        nold = n;
        for (int i = 0; i < nold; i++) {
            told[i] = type[i];
        }
        first = false;
        readcart(ffile,first);
    }

    // perform any final tasks before program exit
    ffile.close();
    if (!test) final();
}

void resizeNumDer() {
    ndesum.allocate(3*n);
    ndeb.allocate(3*n);
    ndea.allocate(3*n);
    ndeba.allocate(3*n);
    ndeub.allocate(3*n);
    ndeaa.allocate(3*n);
    ndeopb.allocate(3*n);
    ndeopd.allocate(3*n);
    ndeid.allocate(3*n);
    ndeit.allocate(3*n);
    ndet.allocate(3*n);
    ndept.allocate(3*n);
    ndebt.allocate(3*n);
    ndeat.allocate(3*n);
    ndett.allocate(3*n);
    ndev.allocate(3*n);
    nder.allocate(3*n);
    ndedsp.allocate(3*n);
    ndec.allocate(3*n);
    ndecd.allocate(3*n);
    nded.allocate(3*n);
    ndem.allocate(3*n);
    ndep.allocate(3*n);
    ndect.allocate(3*n);
    nderxf.allocate(3*n);
    ndes.allocate(3*n);
    ndelf.allocate(3*n);
    ndeg.allocate(3*n);
    ndex.allocate(3*n);
}
}
