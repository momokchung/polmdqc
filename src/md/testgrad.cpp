// Author: Moses KJ Chung
// Year:   2024

#include "atoms.h"
#include "calcMode.h"
#include "deriv.h"
#include "energi.h"
#include "energy.h"
#include "files.h"
#include "getcart.h"
#include "gettext.h"
#include "inform.h"
#include "initial.h"
#include "libfunc.h"
#include "mechanic.h"
#include "nextarg.h"
#include "output.h"
#include "readcart.cpp"
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
            printf("\n Analysis for Archive Structure :        %8d\n", frame);
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
        if (doanalyt) {
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
        if (dofull) {
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
            if (donumer and use[i]) {
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
                    ndesum[i][j] = (f - f0) / eps;
                    ndeb[i][j] = (eb - eb0) / eps;
                    ndea[i][j] = (ea - ea0) / eps;
                    ndeba[i][j] = (eba - eba0) / eps;
                    ndeub[i][j] = (eub - eub0) / eps;
                    ndeaa[i][j] = (eaa - eaa0) / eps;
                    ndeopb[i][j] = (eopb - eopb0) / eps;
                    ndeopd[i][j] = (eopd - eopd0) / eps;
                    ndeid[i][j] = (eid - eid0) / eps;
                    ndeit[i][j] = (eit - eit0) / eps;
                    ndet[i][j] = (et - et0) / eps;
                    ndept[i][j] = (ept - ept0) / eps;
                    ndebt[i][j] = (ebt - ebt0) / eps;
                    ndeat[i][j] = (eat - eat0) / eps;
                    ndett[i][j] = (ett - ett0) / eps;
                    ndev[i][j] = (ev - ev0) / eps;
                    nder[i][j] = (er - er0) / eps;
                    ndedsp[i][j] = (edsp - edsp0) / eps;
                    ndec[i][j] = (ec - ec0) / eps;
                    ndecd[i][j] = (ecd - ecd0) / eps;
                    nded[i][j] = (ed - ed0) / eps;
                    ndem[i][j] = (em - em0) / eps;
                    ndep[i][j] = (ep - ep0) / eps;
                    ndect[i][j] = (ect - ect0) / eps;
                    nderxf[i][j] = (erxf - erxf0) / eps;
                    ndes[i][j] = (es - es0) / eps;
                    ndelf[i][j] = (elf - elf0) / eps;
                    ndeg[i][j] = (eg - eg0) / eps;
                    ndex[i][j] = (ex - ex0) / eps;
                }
            }

            // print analytical gradients of each energy term for each atom
            if (dofull and use[i]) {
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
                                i+1,deb[i][j],dea[i][j],deba[i][j],deub[i][j],
                                axis[j],deaa[i][j],deopb[i][j],deopd[i][j],deid[i][j],
                                deit[i][j],det[i][j],dept[i][j],debt[i][j],
                                deat[i][j],dett[i][j],dev[i][j],der[i][j],
                                dedsp[i][j],dec[i][j],decd[i][j],ded[i][j],
                                dem[i][j],dep[i][j],dect[i][j],derxf[i][j],
                                des[i][j],delf[i][j],deg[i][j],dex[i][j]);
                        }
                        else if (digits >= 6) {
                            printf("\n%6d%14.6f%14.6f%14.6f%14.6f%14.6f"
                               "\n     %c%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n Anlyt%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f\n",
                                i+1,deb[i][j],dea[i][j],deba[i][j],deub[i][j],deaa[i][j],
                                axis[j],deopb[i][j],deopd[i][j],deid[i][j],deit[i][j],det[i][j],
                                dept[i][j],debt[i][j],deat[i][j],dett[i][j],dev[i][j],
                                der[i][j],dedsp[i][j],dec[i][j],decd[i][j],ded[i][j],
                                dem[i][j],dep[i][j],dect[i][j],derxf[i][j],des[i][j],
                                delf[i][j],deg[i][j],dex[i][j]);
                        }
                        else {
                            printf("\n%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                               "\n     %c%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n Anlyt%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f\n",
                                i+1,deb[i][j],dea[i][j],deba[i][j],deub[i][j],deaa[i][j],deopb[i][j],
                                axis[j],deopd[i][j],deid[i][j],deit[i][j],det[i][j],dept[i][j],debt[i][j],
                                deat[i][j],dett[i][j],dev[i][j],der[i][j],dedsp[i][j],dec[i][j],
                                decd[i][j],ded[i][j],dem[i][j],dep[i][j],dect[i][j],derxf[i][j],
                                des[i][j],delf[i][j],deg[i][j],dex[i][j]);
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
                                i+1,ndeb[i][j],ndea[i][j],ndeba[i][j],ndeub[i][j],
                                axis[j],ndeaa[i][j],ndeopb[i][j],ndeopd[i][j],ndeid[i][j],
                                ndeit[i][j],ndet[i][j],ndept[i][j],ndebt[i][j],
                                ndeat[i][j],ndett[i][j],ndev[i][j],nder[i][j],
                                ndedsp[i][j],ndec[i][j],ndecd[i][j],nded[i][j],
                                ndem[i][j],ndep[i][j],ndect[i][j],nderxf[i][j],
                                ndes[i][j],ndelf[i][j],ndeg[i][j],ndex[i][j]);
                        }
                        else if (digits >= 6) {
                            printf("\n%6d%14.6f%14.6f%14.6f%14.6f%14.6f"
                               "\n     %c%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n Numer%14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f%14.6f%14.6f"
                                "\n      %14.6f%14.6f%14.6f\n",
                                i+1,ndeb[i][j],ndea[i][j],ndeba[i][j],ndeub[i][j],ndeaa[i][j],
                                axis[j],ndeopb[i][j],ndeopd[i][j],ndeid[i][j],ndeit[i][j],ndet[i][j],
                                ndept[i][j],ndebt[i][j],ndeat[i][j],ndett[i][j],ndev[i][j],
                                nder[i][j],ndedsp[i][j],ndec[i][j],ndecd[i][j],nded[i][j],
                                ndem[i][j],ndep[i][j],ndect[i][j],nderxf[i][j],ndes[i][j],
                                ndelf[i][j],ndeg[i][j],ndex[i][j]);
                        }
                        else {
                            printf("\n%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                               "\n     %c%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n Numer%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f%12.4f%12.4f"
                                "\n      %12.4f%12.4f%12.4f%12.4f\n",
                                i+1,ndeb[i][j],ndea[i][j],ndeba[i][j],ndeub[i][j],ndeaa[i][j],ndeopb[i][j],
                                axis[j],ndeopd[i][j],ndeid[i][j],ndeit[i][j],ndet[i][j],ndept[i][j],ndebt[i][j],
                                ndeat[i][j],ndett[i][j],ndev[i][j],nder[i][j],ndedsp[i][j],ndec[i][j],
                                ndecd[i][j],nded[i][j],ndem[i][j],ndep[i][j],ndect[i][j],nderxf[i][j],
                                ndes[i][j],ndelf[i][j],ndeg[i][j],ndex[i][j]);
                        }
                    }
                }
            }
        }

        // print the total gradient components for each atom
        if (doanalyt or donumer) {
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
        totnorm = 0.;
        ntotnorm = 0.;
        for (int i = 0; i < n; i++) {
            if (doanalyt and use[i]) {
                denorm[i] = REAL_POW(desum[i][0],2) + REAL_POW(desum[i][1],2) + REAL_POW(desum[i][2],2);
                totnorm = totnorm + denorm[i];
                denorm[i] = REAL_SQRT(denorm[i]);
                if (digits >= 8) {
                    printf(" Anlyt%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,desum[i][0],desum[i][1],desum[i][2],denorm[i]);
                }
                else if (digits >= 6) {
                    printf(" Anlyt  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,desum[i][0],desum[i][1],desum[i][2],denorm[i]);
                }
                else {
                    printf(" Anlyt  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,desum[i][0],desum[i][1],desum[i][2],denorm[i]);
                }
            }
            if (donumer and use[i]) {
                ndenorm[i] = REAL_POW(ndesum[i][0],2) + REAL_POW(ndesum[i][1],2) + REAL_POW(ndesum[i][2],2);
                ntotnorm = ntotnorm + ndenorm[i];
                ndenorm[i] = REAL_SQRT(ndenorm[i]);
                if (digits >= 8) {
                    printf(" Numer%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,ndesum[i][0],ndesum[i][1],ndesum[i][2],ndenorm[i]);
                }
                else if (digits >= 6) {
                    printf(" Numer  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,ndesum[i][0],ndesum[i][1],ndesum[i][2],ndenorm[i]);
                }
                else {
                    printf(" Numer  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,ndesum[i][0],ndesum[i][1],ndesum[i][2],ndenorm[i]);
                }
            }
        }

        // print the total norm for the analytical gradient
        if (doanalyt or donumer) {
            printf("\n Total Gradient Norm and RMS Gradient per Atom :\n");
        }
        if (doanalyt) {
            totnorm = REAL_SQRT(totnorm);
            if (digits >= 8) {
                printf("\n Anlyt      Total Gradient Norm Value      %20.8f\n", totnorm);
            }
            else if (digits >= 6) {
                printf("\n Anlyt      Total Gradient Norm Value      %18.6f\n", totnorm);
            }
            else {
                printf("\n Anlyt      Total Gradient Norm Value      %16.4f\n", totnorm);
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

        // print the rms per atom norm for the analytical gradient
        if (doanalyt) {
            rms = totnorm / REAL_SQRT(static_cast<real>(nuse));
            if (digits >= 8) {
                printf("\n Anlyt      RMS Gradient over All Atoms    %20.8f\n", rms);
            }
            else if (digits >= 6) {
                printf("\n Anlyt      RMS Gradient over All Atoms    %18.6f\n", rms);
            }
            else {
                printf("\n Anlyt      RMS Gradient over All Atoms    %16.4f\n", rms);
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
}

void resizeNumDer() {
    ndesum.resize(n, std::vector<real>(3));
    ndeb.resize(n, std::vector<real>(3));
    ndea.resize(n, std::vector<real>(3));
    ndeba.resize(n, std::vector<real>(3));
    ndeub.resize(n, std::vector<real>(3));
    ndeaa.resize(n, std::vector<real>(3));
    ndeopb.resize(n, std::vector<real>(3));
    ndeopd.resize(n, std::vector<real>(3));
    ndeid.resize(n, std::vector<real>(3));
    ndeit.resize(n, std::vector<real>(3));
    ndet.resize(n, std::vector<real>(3));
    ndept.resize(n, std::vector<real>(3));
    ndebt.resize(n, std::vector<real>(3));
    ndeat.resize(n, std::vector<real>(3));
    ndett.resize(n, std::vector<real>(3));
    ndev.resize(n, std::vector<real>(3));
    nder.resize(n, std::vector<real>(3));
    ndedsp.resize(n, std::vector<real>(3));
    ndec.resize(n, std::vector<real>(3));
    ndecd.resize(n, std::vector<real>(3));
    nded.resize(n, std::vector<real>(3));
    ndem.resize(n, std::vector<real>(3));
    ndep.resize(n, std::vector<real>(3));
    ndect.resize(n, std::vector<real>(3));
    nderxf.resize(n, std::vector<real>(3));
    ndes.resize(n, std::vector<real>(3));
    ndelf.resize(n, std::vector<real>(3));
    ndeg.resize(n, std::vector<real>(3));
    ndex.resize(n, std::vector<real>(3));
}
}
