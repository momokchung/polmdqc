////////////////////////////////////////////////////////
//                                                    //
//  readprm.cpp  --  input of force field parameters  //
//                                                    //
////////////////////////////////////////////////////////

// "readprm" processes the potential energy parameter file
// in order to define the default force field parameters


#include "fatal.h"
#include "fields.h"
#include "getnumb.h"
#include "getstring.h"
#include "gettext.h"
#include "getword.h"
#include "kanang.h"
#include "kangs.h"
#include "kantor.h"
#include "katoms.h"
#include "kbonds.h"
#include "kcflux.h"
#include "kchrge.h"
#include "kcpen.h"
#include "kctrn.h"
#include "kdipol.h"
#include "kdsp.h"
#include "kexpl.h"
#include "khbond.h"
#include "kiprop.h"
#include "kitors.h"
#include "kmulti.h"
#include "kopbnd.h"
#include "kopdst.h"
#include "korbs.h"
#include "kpitor.h"
#include "kpolpr.h"
#include "kpolr.h"
#include "krepl.h"
#include "ksolut.h"
#include "kstbnd.h"
#include "ksttor.h"
#include "ktorsn.h"
#include "ktrtor.h"
#include "kurybr.h"
#include "kvdws.h"
#include "kvdwpr.h"
#include "merck.h"
#include "numeral.h"
#include "params.h"
#include "prmkey.h"
#include "readprm.h"
#include "solute.h"
#include "sort.h"
#include "torphase.h"
#include "upcase.h"
#include <algorithm>
#include <sstream>

void readprm()
{
    int iprm;
    int ia,ib,ic,id;
    int ie,ig,ih,ii,ij;
    int imin;
    int size,next;
    int cls,atn,lig;
    int nx,ny,nxy;
    int bt,at,sbt,tt;
    int ilpr;
    int ft[6],pg[maxval];
    double wght,rd;
    double ep,rdn;
    double spr,apr,epr;
    double cdp,adp;
    double an1,an2,an3;
    double ba1,ba2;
    double aa1,aa2,aa3;
    double bt1,bt2,bt3;
    double bt4,bt5,bt6;
    double bt7,bt8,bt9;
    double at1,at2,at3;
    double at4,at5,at6;
    double an,pr,ds,dk;
    double vd,cg,dp,ps;
    double fc,bd,dl;
    double pt,pel,pal;
    double pol,thl,thd;
    double kpr,ppr,dpr;
    double ctrn,atrn;
    double cfb,cfb1,cfb2;
    double cfa1,cfa2;
    double pbrd,csrd,gkrd;
    double el,iz,rp;
    double ss,ts;
    double abc,cba;
    double gi,alphi;
    double nni,factor;
    double vt[6],st[6];
    double pl[13];
    bool header,swap;
    char da1;
    std::string pa,pb,pc;
    std::string pd,pe;
    std::string axt;
    std::string keyword;
    std::string text;
    std::string record;
    std::string string;
    std::istringstream iss;

    // initialize the counters for some parameter types
    int nb = 0;
    int nb5 = 0;
    int nb4 = 0;
    int nb3 = 0;
    int nel = 0;
    int na = 0;
    int na5 = 0;
    int na4 = 0;
    int na3 = 0;
    int nap = 0;
    int naf = 0;
    int nsb = 0;
    int nu = 0;
    int nopb = 0;
    int nopd = 0;
    int ndi = 0;
    int nti = 0;
    int nt = 0;
    int nt5 = 0;
    int nt4 = 0;
    int npt = 0;
    int nbt = 0;
    int nat = 0;
    int ntt = 0;
    int nvp = 0;
    int nhb = 0;
    int nd = 0;
    int nd5 = 0;
    int nd4 = 0;
    int nd3 = 0;
    int nmp = 0;
    int npp = 0;
    int ncfb = 0;
    int ncfa = 0;
    int npi = 0;
    int npi5 = 0;
    int npi4 = 0;

    // number of characters in an atom number text string
    size = 4;

    // set blank line header before echoed comment lines
    header = true;

    // process each line of the parameter file, first
    // extract the keyword at the start of each line
    iprm = 0;
    while (iprm < nprm) {
        record = prmline[iprm];
        next = 0;
        gettext(record,keyword,next);
        upcase(keyword);
        string = record.substr(next);
        iss.clear();
        iss.str(string);

        // check for a force field modification keyword
        prmkey(record);

        // comment line to be echoed to the output
        if (keyword == "ECHO") {
            if (header) {
               header = false;
               printf("\n");
            }
            printf("%s\n", string.c_str());
        }

        // atom type definitions and parameters
        else if (keyword == "ATOM") {
            ia = -1;
            cls = -1;
            atn = 0;
            wght = 0.;
            lig = 0;
            getnumb(record,ia,next);
            getnumb(record,cls,next);
            ia--;
            cls--;
            if (cls == -1)  cls = ia;
            atmcls[ia] = cls;
            if (ia >= maxtyp) {
                printf("\n READPRM  --  Too many Atom Types; Increase MAXTYP\n");
                fatal();
            }
            else if (cls >= maxclass) {
                printf("\n READPRM  --  Too many Atom Classes; Increase MAXCLASS\n");
                fatal();
            }
            if (ia != -1) {
                gettext(record,symbol[ia],next);
                getstring(record,describe[ia],next);
                string = record.substr(next);
                iss.clear();
                iss.str(string);
                iss >> atn >> wght >> lig;
                atmnum[ia] = atn;
                weight[ia] = wght;
                ligand[ia] = lig;
            }
        }

        // bond stretching parameters
        else if (keyword == "BOND") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> ib >> fc >> bd;
            ia--;
            ib--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            if (ia <= ib) {
               kb[nb] = pa + pb;
            }
            else {
               kb[nb] = pb + pa;
            }
            bcon[nb] = fc;
            blen[nb] = bd;
            nb++;
        }

        // bond stretching parameters for 5-membered rings
        else if (keyword == "BOND5") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> ib >> fc >> bd;
            ia--;
            ib--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            if (ia <= ib) {
               kb5[nb5] = pa + pb;
            }
            else {
               kb5[nb5] = pb + pa;
            }
            bcon5[nb5] = fc;
            blen5[nb5] = bd;
            nb5++;
        }

        // bond stretching parameters for 4-membered rings
        else if (keyword == "BOND4") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> ib >> fc >> bd;
            ia--;
            ib--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            if (ia <= ib) {
               kb4[nb4] = pa + pb;
            }
            else {
               kb4[nb4] = pb + pa;
            }
            bcon4[nb4] = fc;
            blen4[nb4] = bd;
            nb4++;
        }

        // bond stretching parameters for 3-membered rings
        else if (keyword == "BOND3") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> ib >> fc >> bd;
            ia--;
            ib--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            if (ia <= ib) {
               kb3[nb3] = pa + pb;
            }
            else {
               kb3[nb3] = pb + pa;
            }
            bcon3[nb3] = fc;
            blen3[nb3] = bd;
            nb3++;
        }

        // electronegativity bond length correction parameters
        else if (keyword == "ELECTNEG") {
            ia = -1;
            ib = -1;
            ic = 0;
            dl = 0.;
            iss >> ia >> ib >> ic >> dl;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               kel[nel] = pa+pb+pc;
            }
            else {
               kel[nel] = pc+pb+pa;
            }
            dlen[nel] = dl;
            nel++;
        }

        // bond angle bending parameters
        else if (keyword == "ANGLE") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            an3 = 0.;
            iss >> ia >> ib >> ic >> fc >> an1 >> an2 >> an3;
            ia--;
            ib--;
            ic--;
            pa = numeral (ia, size);
            pb = numeral (ib, size);
            pc = numeral (ic, size);
            if (ia <= ic) {
               ka[na] = pa + pb + pc;
            }
            else {
               ka[na] = pc + pb + pa;
            }
            acon[na] = fc;
            ang[na][0] = an1;
            ang[na][1] = an2;
            ang[na][2] = an3;
            na++;
        }

        // angle bending parameters for 5-membered rings
        else if (keyword == "ANGLE5") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            an3 = 0.;
            iss >> ia >> ib >> ic >> fc >> an1 >> an2 >> an3;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               ka5[na5] = pa + pb + pc;
            }
            else {
               ka5[na5] = pc + pb + pa;
            }
            acon5[na5] = fc;
            ang5[na5][0] = an1;
            ang5[na5][1] = an2;
            ang5[na5][2] = an3;
            na5++;
        }

        // angle bending parameters for 4-membered rings
        else if (keyword == "ANGLE4") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            an3 = 0.;
            iss >> ia >> ib >> ic >> fc >> an1 >> an2 >> an3;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               ka4[na4] = pa + pb + pc;
            }
            else {
               ka4[na4] = pc + pb + pa;
            }
            acon4[na4] = fc;
            ang4[na4][0] = an1;
            ang4[na4][1] = an2;
            ang4[na4][2] = an3;
            na4++;
        }

        // angle bending parameters for 3-membered rings
        else if (keyword == "ANGLE3") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            an3 = 0.;
            iss >> ia >> ib >> ic >> fc >> an1 >> an2 >> an3;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               ka3[na3] = pa + pb + pc;
            }
            else {
               ka3[na3] = pc + pb + pa;
            }
            acon3[na3] = fc;
            ang3[na3][0] = an1;
            ang3[na3][1] = an2;
            ang3[na3][2] = an3;
            na3++;
        }

        // in-plane projected angle bending parameters
        else if (keyword == "ANGLEP") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an1 = 0.;
            an2 = 0.;
            iss >> ia >> ib >> ic >> fc >> an1 >> an2;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               kap[nap] = pa + pb + pc;
            }
            else {
               kap[nap] = pc + pb + pa;
            }
            aconp[nap] = fc;
            angp[nap][0] = an1;
            angp[nap][1] = an2;
            nap++;
        }

        // Fourier bond angle bending parameters
        else if (keyword == "ANGLEF") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            an = 0.;
            pr = 0.;
            iss >> ia >> ib >> ic >> fc >> an >> pr;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               kaf[naf] = pa + pb + pc;
            }
            else {
               kaf[naf] = pc + pb + pa;
            }
            aconf[naf] = fc;
            angf[naf][0] = an;
            angf[naf][1] = pr;
            naf++;
        }

        // stretch-bend parameters
        else if (keyword == "STRBND") {
            ia = -1;
            ib = -1;
            ic = -1;
            ba1 = 0.;
            ba2 = 0.;
            iss >> ia >> ib >> ic >> ba1 >> ba2;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               ksb[nsb] = pa + pb + pc;
               stbn[nsb][0] = ba1;
               stbn[nsb][1] = ba2;
            }
            else {
               ksb[nsb] = pc + pb + pa;
               stbn[nsb][0] = ba2;
               stbn[nsb][1] = ba1;
            }
            nsb++;
        }

        // Urey-Bradley parameters
        else if (keyword == "UREYBRAD") {
            ia = -1;
            ib = -1;
            ic = -1;
            fc = 0.;
            ds = 0.;
            iss >> ia >> ib >> ic >> fc >> ds;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            if (ia <= ic) {
               ku[nu] = pa + pb + pc;
            }
            else {
               ku[nu] = pc + pb + pa;
            }
            ucon[nu] = fc;
            dst13[nu] = ds;
            nu++;
        }

        // angle-angle parameters
        else if (keyword == "ANGANG") {
            ia = -1;
            aa1 = 0.;
            aa2 = 0.;
            aa3 = 0.;
            iss >> ia >> aa1 >> aa2 >> aa3;
            ia--;
            if (ia != -1) {
               anan[ia][0] = aa1;
               anan[ia][1] = aa2;
               anan[ia][2] = aa3;
            }
        }

        // out-of-plane bend parameters
        else if (keyword == "OPBEND") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            fc = 0.;
            iss >> ia >> ib >> ic >> id >> fc;
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ic <= id) {
               kopb[nopb] = pa + pb + pc + pd;
            }
            else {
               kopb[nopb] = pa + pb + pd + pc;
            }
            opbn[nopb] = fc;
            nopb++;
        }

        // out-of-plane distance parameters
        else if (keyword == "OPDIST") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            fc = 0.;
            iss >> ia >> ib >> ic >> id >> fc;
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            imin = std::min({ib,ic,id});
            if (ib == imin) {
                if (ic <= id) {
                  kopd[nopd] = pa + pb + pc + pd;
                }
                else {
                  kopd[nopd] = pa + pb + pd + pc;
                }
            }
            else if (ic == imin) {
                if (ib <= id) {
                  kopd[nopd] = pa + pc + pb + pd;
                }
                else {
                  kopd[nopd] = pa + pc + pd + pb;
                }
            }
            else if (id == imin) {
                if (ib <= ic) {
                  kopd[nopd] = pa + pd + pb + pc;
                }
                else {
                  kopd[nopd] = pa + pd + pc + pb;
                }
            }
            opds[nopd] = fc;
            nopd++;
        }

        // improper dihedral parameters
        else if (keyword == "IMPROPER") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            dk = 0.;
            vd = 0.;
            iss >> ia >> ib >> ic >> id >> dk >> vd;
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            kdi[ndi] = pa + pb + pc + pd;
            dcon[ndi] = dk;
            tdi[ndi] = vd;
            ndi++;
        }

        // improper torsional parameters
        else if (keyword == "IMPTORS") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            for (int i = 0; i < 6; i++) {
               vt[i] = 0.;
               st[i] = 0.;
               ft[i] = -1;
            }
            iss >> ia >> ib >> ic >> id;
            ia--;
            ib--;
            ic--;
            id--;
            for (int i = 0; i < 6; i++) {
                iss >> vt[i] >> st[i] >> ft[i];
                ft[i]--;
            }
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            kti[nti] = pa + pb + pc + pd;
            torphase(ft,vt,st);
            ti1[nti][0] = vt[0];
            ti1[nti][1] = st[0];
            ti2[nti][0] = vt[1];
            ti2[nti][1] = st[1];
            ti3[nti][0] = vt[2];
            ti3[nti][1] = st[2];
            nti++;
        }

        // torsional parameters
        else if (keyword == "TORSION") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            for (int i = 0; i < 6; i++) {
               vt[i] = 0.;
               st[i] = 0.;
               ft[i] = -1;
            }
            iss >> ia >> ib >> ic >> id;
            ia--;
            ib--;
            ic--;
            id--;
            for (int i = 0; i < 6; i++) {
                iss >> vt[i] >> st[i] >> ft[i];
                ft[i]--;
            }
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ib < ic) {
                kt[nt] = pa + pb + pc + pd;
            }
            else if (ic < ib) {
                kt[nt] = pd + pc + pb + pa;
            }
            else if (ia <= id) {
                kt[nt] = pa + pb + pc + pd;
            }
            else if (id < ia) {
                kt[nt] = pd + pc + pb + pa;
            }
            torphase(ft,vt,st);
            t1[nt][0] = vt[0];
            t1[nt][1] = st[0];
            t2[nt][0] = vt[1];
            t2[nt][1] = st[1];
            t3[nt][0] = vt[2];
            t3[nt][1] = st[2];
            t4[nt][0] = vt[3];
            t4[nt][1] = st[3];
            t5[nt][0] = vt[4];
            t5[nt][1] = st[4];
            t6[nt][0] = vt[5];
            t6[nt][1] = st[5];
            nt++;
        }

        // torsional parameters for 5-membered rings
        else if (keyword == "TORSION5") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            for (int i = 0; i < 6; i++) {
               vt[i] = 0.;
               st[i] = 0.;
               ft[i] = -1;
            }
            iss >> ia >> ib >> ic >> id;
            ia--;
            ib--;
            ic--;
            id--;
            for (int i = 0; i < 6; i++) {
                iss >> vt[i] >> st[i] >> ft[i];
                ft[i]--;
            }
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ib < ic) {
                kt5[nt5] = pa + pb + pc + pd;
            }
            else if (ic < ib) {
                kt5[nt5] = pd + pc + pb + pa;
            }
            else if (ia <= id) {
                kt5[nt5] = pa + pb + pc + pd;
            }
            else if (id < ia) {
                kt5[nt5] = pd + pc + pb + pa;
            }
            torphase(ft,vt,st);
            t15[nt5][0] = vt[0];
            t15[nt5][1] = st[0];
            t25[nt5][0] = vt[1];
            t25[nt5][1] = st[1];
            t35[nt5][0] = vt[2];
            t35[nt5][1] = st[2];
            t45[nt5][0] = vt[3];
            t45[nt5][1] = st[3];
            t55[nt5][0] = vt[4];
            t55[nt5][1] = st[4];
            t65[nt5][0] = vt[5];
            t65[nt5][1] = st[5];
            nt5++;
        }

        // torsional parameters for 4-membered rings
        else if (keyword == "TORSION4") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            for (int i = 0; i < 6; i++) {
               vt[i] = 0.;
               st[i] = 0.;
               ft[i] = -1;
            }
            iss >> ia >> ib >> ic >> id;
            ia--;
            ib--;
            ic--;
            id--;
            for (int i = 0; i < 6; i++) {
                iss >> vt[i] >> st[i] >> ft[i];
                ft[i]--;
            }
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ib < ic) {
                kt4[nt4] = pa + pb + pc + pd;
            }
            else if (ic < ib) {
                kt4[nt4] = pd + pc + pb + pa;
            }
            else if (ia <= id) {
                kt4[nt4] = pa + pb + pc + pd;
            }
            else if (id < ia) {
                kt4[nt4] = pd + pc + pb + pa;
            }
            torphase(ft,vt,st);
            t14[nt4][0] = vt[0];
            t14[nt4][1] = st[0];
            t24[nt4][0] = vt[1];
            t24[nt4][1] = st[1];
            t34[nt4][0] = vt[2];
            t34[nt4][1] = st[2];
            t44[nt4][0] = vt[3];
            t44[nt4][1] = st[3];
            t54[nt4][0] = vt[4];
            t54[nt4][1] = st[4];
            t64[nt4][0] = vt[5];
            t64[nt4][1] = st[5];
            nt4++;
        }

        // pi-system torsion parameters
        else if (keyword == "PITORS") {
            ia = -1;
            ib = -1;
            pt = 0.;
            iss >> ia >> ib >> pt;
            ia--;
            ib--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            if (ia <= ib) {
               kpt[npt] = pa + pb;
            }
            else {
               kpt[npt] = pb + pa;
            }
            ptcon[npt] = pt;
            npt++;
        }

        // stretch-torsion parameters
        else if (keyword == "STRTORS") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            bt1 = 0.;
            bt2 = 0.;
            bt3 = 0.;
            bt4 = 0.;
            bt5 = 0.;
            bt6 = 0.;
            bt7 = 0.;
            bt8 = 0.;
            bt9 = 0.;
            iss >> ia >> ib >> ic >> id >> bt1 >> bt2 >> bt3 >> bt4 >> bt5 >> bt6 >> bt7 >> bt8 >> bt9;
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ib < ic) {
                kbt[nbt] = pa + pb + pc + pd;
                swap = false;
            }
            else if (ic < ib) {
                kbt[nbt] = pd + pc + pb + pa;
                swap = true;
            }
            else if (ia <= id) {
                kbt[nbt] = pa + pb + pc + pd;
                swap = false;
            }
            else if (id < ia) {
                kbt[nbt] = pd + pc + pb + pa;
                swap = true;
            }
            btcon[nbt][3] = bt4;
            btcon[nbt][4] = bt5;
            btcon[nbt][5] = bt6;
            if (swap) {
                btcon[nbt][0] = bt7;
                btcon[nbt][1] = bt8;
                btcon[nbt][2] = bt9;
                btcon[nbt][6] = bt1;
                btcon[nbt][7] = bt2;
                btcon[nbt][8] = bt3;
            }
            else {
                btcon[nbt][0] = bt1;
                btcon[nbt][1] = bt2;
                btcon[nbt][2] = bt3;
                btcon[nbt][6] = bt7;
                btcon[nbt][7] = bt8;
                btcon[nbt][8] = bt9;
            }
            nbt++;
        }

        // angle-torsion parameters
        else if (keyword == "ANGTORS") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            at1 = 0.;
            at2 = 0.;
            at3 = 0.;
            at4 = 0.;
            at5 = 0.;
            at6 = 0.;
            iss >> ia >> ib >> ic >> id >> at1 >> at2 >> at3 >> at4 >> at5 >> at6;
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            if (ib < ic) {
               kat[nat] = pa + pb + pc + pd;
               swap = false;
            }
            else if (ic < ib) {
               kat[nat] = pd + pc + pb + pa;
               swap = true;
            }
            else if (ia <= id) {
               kat[nat] = pa + pb + pc + pd;
               swap = false;
            }
            else if (id < ia) {
               kat[nat] = pd + pc + pb + pa;
               swap = true;
            }
            if (swap) {
               atcon[nat][0] = at4;
               atcon[nat][1] = at5;
               atcon[nat][2] = at6;
               atcon[nat][3] = at1;
               atcon[nat][4] = at2;
               atcon[nat][5] = at3;
            }
            else {
               atcon[nat][0] = at1;
               atcon[nat][1] = at2;
               atcon[nat][2] = at3;
               atcon[nat][3] = at4;
               atcon[nat][4] = at5;
               atcon[nat][5] = at6;
            }
            nat++;
        }

        // torsion-torsion parameters
        else if (keyword == "TORTORS") {
            ia = -1;
            ib = -1;
            ic = -1;
            id = -1;
            ie = -1;
            nx = 0;
            ny = 0;
            nxy = 0;
            iss >> ia >> ib >> ic >> id >> ie >> nx >> ny;
            ia--;
            ib--;
            ic--;
            id--;
            ie--;
            nxy = nx * ny;
            std::vector<double> tx(nxy, 0.);
            std::vector<double> ty(nxy, 0.);
            std::vector<double> tf(nxy, 0.);
            for (int i = 0; i < nxy; i++) {
                iprm++;
                record = prmline[iprm];
                iss.clear();
                iss.str(record);
                iss >> tx[i] >> ty[i] >> tf[i];
            }
            pa = numeral(ia, size);
            pb = numeral(ib, size);
            pc = numeral(ic, size);
            pd = numeral(id, size);
            pe = numeral(ie, size);
            ktt[ntt] = pa + pb + pc + pd + pe;
            nx = nxy;
            sort(tx);
            ny = nxy;
            sort(ty);
            tnx[ntt] = nx;
            tny[ntt] = ny;
            for (int i = 0; i < nx; i++) {
                ttx[ntt][i] = tx[i];
            }
            for (int i = 0; i < ny; i++) {
                tty[ntt][i] = ty[i];
            }
            for (int i = 0; i < nxy; i++) {
                tbf[ntt][i] = tf[i];
            }
            ntt++;
        }

        // van der Waals parameters for individual atom types
        else if (keyword == "VDW") {
            ia = -1;
            rd = 0.;
            ep = 0.;
            rdn = 0.;
            iss >> ia >> rd >> ep >> rdn;
            ia--;
            if (ia != -1) {
                rad[ia] = rd;
                eps[ia] = ep;
                reduct[ia] = rdn;
            }
        }

        // van der Waals 1-4 parameters for individual atom types
        else if (keyword == "VDW14") {
            ia = -1;
            rd = 0.;
            ep = 0.;
            iss >> ia >> rd >> ep;
            ia--;
            if (ia != -1) {
               rad4[ia] = rd;
               eps4[ia] = ep;
            }
        }

        // van der Waals parameters for specific atom pairs
        else if (keyword == "VDWPAIR" or keyword == "VDWPR") {
            ia = -1;
            ib = -1;
            rd = 0.;
            ep = 0.;
            iss >> ia >> ib >> rd >> ep;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kvpr[nvp] = pa + pb;
            }
            else {
               kvpr[nvp] = pb + pa;
            }
            radpr[nvp] = rd;
            epspr[nvp] = ep;
            nvp++;
        }

        // van der Waals parameters for hydrogen bonding pairs
        else if (keyword == "HBOND") {
            ia = -1;
            ib = -1;
            rd = 0.;
            ep = 0.;
            iss >> ia >> ib >> rd >> ep;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               khb[nhb] = pa + pb;
            }
            else {
               khb[nhb] = pb + pa;
            }
            radhb[nhb] = rd;
            epshb[nhb] = ep;
            nhb++;
        }

        // Pauli repulsion parameters
        else if (keyword == "REPULSION") {
            ia = -1;
            spr = 0.;
            apr = 0.;
            epr = 0.;
            iss >> ia >> spr >> apr >> epr;
            ia--;
            if (ia != -1) {
               prsiz[ia] = spr;
               prdmp[ia] = apr;
               prele[ia] = -std::abs(epr);
            }
        }

        // damped dispersion parameters
        else if (keyword == "DISPERSION") {
            ia = -1;
            cdp = 0.;
            adp = 0.;
            iss >> ia >> cdp >> adp;
            ia--;
            if (ia != -1) {
               dspsix[ia] = cdp;
               dspdmp[ia] = adp;
            }
        }

        // atomic partial charge parameters
        else if (keyword == "CHARGE") {
            ia = -1;
            cg = 0.;
            iss >> ia >> cg;
            ia--;
            if (ia != -1)  chg[ia] = cg;
        }

        // bond dipole moment parameters
        else if (keyword == "DIPOLE") {
            ia = -1;
            ib = -1;
            dp = 0.0;
            ps = 0.5;
            iss >> ia >> ib >> dp >> ps;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kd[nd] = pa + pb;
               dpl[nd] = dp;
               pos[nd] = ps;
            }
            else {
               kd[nd] = pb + pa;
               dpl[nd] = -dp;
               pos[nd] = 1. - ps;
            }
            nd++;
        }

        // bond dipole moment parameters for 5-membered rings
        else if (keyword == "DIPOLE5") {
            ia = -1;
            ib = -1;
            dp = 0.0;
            ps = 0.5;
            iss >> ia >> ib >> dp >> ps;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kd5[nd5] = pa + pb;
               dpl5[nd5] = dp;
               pos5[nd5] = ps;
            }
            else {
               kd5[nd5] = pb + pa;
               dpl5[nd5] = -dp;
               pos5[nd5] = 1. - ps;
            }
            nd5++;
        }

        // bond dipole moment parameters for 4-membered rings
        else if (keyword == "DIPOLE4") {
            ia = -1;
            ib = -1;
            dp = 0.0;
            ps = 0.5;
            iss >> ia >> ib >> dp >> ps;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kd4[nd4] = pa + pb;
               dpl4[nd4] = dp;
               pos4[nd4] = ps;
            }
            else {
               kd4[nd4] = pb + pa;
               dpl4[nd4] = -dp;
               pos4[nd4] = 1. - ps;
            }
            nd4++;
        }

        // bond dipole moment parameters for 3-membered rings
        else if (keyword == "DIPOLE3") {
            ia = -1;
            ib = -1;
            dp = 0.0;
            ps = 0.5;
            iss >> ia >> ib >> dp >> ps;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kd3[nd3] = pa + pb;
               dpl3[nd3] = dp;
               pos3[nd3] = ps;
            }
            else {
               kd3[nd3] = pb + pa;
               dpl3[nd3] = -dp;
               pos3[nd3] = 1. - ps;
            }
            nd3++;
        }

        // atomic multipole moment parameters
        else if (keyword == "MULTIPOLE") {
            ia = 0;
            ib = 0;
            ic = 0;
            id = 0;
            axt = "Z-then-X";
            for (int i = 0; i < 13; i++) {
               pl[i] = 0.;
            }
            int numObjects = 0;
            int value;
            while (iss >> value) {
                numObjects++;
            }
            iss.clear();
            iss.seekg(0);
            if (numObjects == 5) {
                iss >> ia >> ib >> ic >> id >> pl[0];
            }
            else if (numObjects == 4) {
                iss >> ia >> ib >> ic >> pl[0];
            }
            else if (numObjects == 3) {
                iss >> ia >> ib >> pl[0];
            }
            else if (numObjects == 2) {
                iss >> ia >> pl[0];
            }
            iprm++;
            record = prmline[iprm];
            iss.clear();
            iss.str(record);
            if (!(iss >> pl[1] >> pl[2] >> pl[3])) goto label_470;
            iprm++;
            record = prmline[iprm];
            iss.clear();
            iss.str(record);
            if (!(iss >> pl[4])) goto label_470;
            iprm++;
            record = prmline[iprm];
            iss.clear();
            iss.str(record);
            if (!(iss >> pl[7] >> pl[8])) goto label_470;
            iprm++;
            record = prmline[iprm];
            iss.clear();
            iss.str(record);
            if (!(iss >> pl[10] >> pl[11] >> pl[12])) goto label_470;
            label_470:
            if (ib == 0) axt = "None";
            if (ib!=0 and ic==0) axt = "Z-Only";
            if (ib<0 or ic<0) axt = "Bisector";
            if (ic<0 and id<0) axt = "Z-Bisect";
            if (std::max({ib,ic,id}) < 0) axt = "3-Fold";
            ib = std::abs(ib);
            ic = std::abs(ic);
            id = std::abs(id);
            ia--;
            ib--;
            ic--;
            id--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            pd = numeral(id,size);
            kmp[nmp] = pa + pb + pc + pd;
            mpaxis[nmp] = axt;
            multip[nmp][0] = pl[0];
            multip[nmp][1] = pl[1];
            multip[nmp][2] = pl[2];
            multip[nmp][3] = pl[3];
            multip[nmp][4] = pl[4];
            multip[nmp][5] = pl[7];
            multip[nmp][6] = pl[10];
            multip[nmp][7] = pl[7];
            multip[nmp][8] = pl[8];
            multip[nmp][9] = pl[11];
            multip[nmp][10] = pl[10];
            multip[nmp][11] = pl[11];
            multip[nmp][12] = pl[12];
            nmp++;
        }

        // charge penetration parameters
        else if (keyword == "CHGPEN") {
            ia = -1;
            pel = 0.;
            pal = 0.;
            iss >> ia >> pel >> pal;
            ia--;
            if (ia != -1) {
                cpele[ia] = std::abs(pel);
                cpalp[ia] = pal;
            }
        }

        // atomic dipole polarizability parameters
        else if (keyword == "POLARIZE") {
            ia = -1;
            pol = 0.;
            thl = 0.;
            thd = 0.;
            for (int i = 0; i < maxval; i++){
               pg[i] = 0;
            }
            string = record;
            getnumb(string,ia,next);
            ia--;
            gettext(string,text,next);
            iss.clear();
            iss.str(text);
            iss >> pol;
            gettext(string,text,next);
            int i = 0;
            getnumb(text,pg[0],i);
            if (pg[0] == 0) {
                iss.clear();
                iss.str(text);
                iss >> thl;
                gettext(string,text,next);
                int i = 0;
                getnumb(text,pg[0],i);
                string = string.substr(next);
                if (pg[0] == 0) {
                    iss.clear();
                    iss.str(text);
                    iss >> thd;
                    iss.clear();
                    iss.str(string);
                    int i = 0;
                    while ((iss >> pg[i]) and i<maxval) {
                        i++;
                    }
                }
                else {
                    iss.clear();
                    iss.str(string);
                    int i = 1;
                    while ((iss >> pg[i]) and i<maxval) {
                        i++;
                    }
                }
            }
            else {
                string = string.substr(next);
                iss.clear();
                iss.str(string);
                int i = 1;
                while ((iss >> pg[i]) and i<maxval) {
                    i++;
                }
            }
            if (ia != -1) {
                polr[ia] = pol;
                athl[ia] = thl;
                dthl[ia] = thd;
                for (int i = 0; i < maxval; i++) {
                    pgrp[ia][i] = pg[i] - 1;
                }
            }
        }

        // polarization parameters for specific atom pairs
        else if (keyword == "POLPAIR") {
            ia = -1;
            ib = -1;
            thl = 0.;
            thd = 0.;
            iss >> ia >> ib >> thl >> thd;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kppr[npp] = pa + pb;
            }
            else {
               kppr[npp] = pb + pa;
            }
            thlpr[npp] = thl;
            thdpr[npp] = thd;
            npp++;
        }

        // exchange polarization parameters
        else if (keyword == "EXCHPOL") {
            ia = -1;
            kpr = 0.;
            ppr = 0.;
            dpr = 0.;
            ilpr = 0;
            iss >> ia >> kpr >> ppr >> dpr >> ilpr;
            ia--;
            if (ia != -1) {
                pepk[ia] = kpr;
                peppre[ia] = ppr;
                pepdmp[ia] = dpr;
                if (ilpr != 0) {
                    pepl[ia] = true;
                }
                else {
                    pepl[ia] = false;
                }
            }
        }

        // charge transfer parameters
        else if (keyword == "CHGTRN") {
            ia = -1;
            ctrn = 0.;
            atrn = 0.;
            iss >> ia >> ctrn >> atrn;
            ia--;
            if (ia != -1) {
                ctchg[ia] = ctrn;
                ctdmp[ia] = atrn;
            }
        }

        // bond charge flux parameters
        else if (keyword == "BNDCFLUX") {
            ia = -1;
            ib = -1;
            cfb = 0.;
            iss >> ia >> ib >> cfb;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia < ib) {
               kcfb[ncfb] = pa + pb;
               cflb[ncfb] = cfb;
            }
            else if (ib < ia) {
               kcfb[ncfb] = pb + pa;
               cflb[ncfb] = -cfb;
            }
            else {
               kcfb[ncfb] = pa + pb;
               cflb[ncfb] = 0.;
            }
            ncfb++;
        }

        // angle charge flux parameters
        else if (keyword == "ANGCFLUX") {
            ia = -1;
            ib = -1;
            ic = -1;
            cfa1 = 0.;
            cfa2 = 0.;
            cfb1 = 0.;
            cfb2 = 0.;
            iss >> ia >> ib >> ic >> cfa1 >> cfa2 >> cfb1 >> cfb2;
            ia--;
            ib--;
            ic--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            pc = numeral(ic,size);
            if (ia <= ic) {
                kcfa[ncfa] = pa + pb + pc;
                cfla[ncfa][0] = cfa1;
                cfla[ncfa][1] = cfa2;
                cflab[ncfa][0] = cfb1;
                cflab[ncfa][1] = cfb2;
            }
            else {
                kcfa[ncfa] = pc + pb + pa;
                cfla[ncfa][0] = cfa2;
                cfla[ncfa][1] = cfa1;
                cflab[ncfa][0] = cfb2;
                cflab[ncfa][1] = cfb1;
            }
            ncfa++;
        }

        // implicit solvation parameters
        else if (keyword == "SOLUTE") {
            ia = -1;
            pbrd = 0.;
            csrd = 0.;
            gkrd = 0.;
            iss >> ia >> pbrd >> csrd >> gkrd;
            ia--;
            if (ia != -1) {
               pbr[ia] = 0.5 * pbrd;
               csr[ia] = 0.5 * csrd;
               gkr[ia] = 0.5 * gkrd;
            }
        }

        // conjugated pisystem atom parameters
        else if (keyword == "PIATOM") {
            ia = -1;
            el = 0.;
            iz = 0.;
            rp = 0.;
            iss >> ia >> el >> iz >> rp;
            ia--;
            if (ia != -1) {
               electron[ia] = el;
               ionize[ia] = iz;
               repulse[ia] = rp;
            }
        }

        // conjugated pisystem bond parameters
        else if (keyword == "PIBOND") {
            ia = -1;
            ib = -1;
            ss = 0.;
            ts = 0.;
            iss >> ia >> ib >> ss >> ts;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kpi[npi] = pa + pb;
            }
            else {
               kpi[npi] = pb + pa;
            }
            sslope[npi] = ss;
            tslope[npi] = ts;
            npi++;
        }

        // conjugated pisystem bond parameters for 5-membered rings
        else if (keyword == "PIBOND5") {
            ia = -1;
            ib = -1;
            ss = 0.;
            ts = 0.;
            iss >> ia >> ib >> ss >> ts;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kpi5[npi5] = pa + pb;
            }
            else {
               kpi5[npi5] = pb + pa;
            }
            sslope5[npi5] = ss;
            tslope5[npi5] = ts;
            npi5++;
        }

        // conjugated pisystem bond parameters for 4-membered rings
        else if (keyword == "PIBOND4") {
            ia = -1;
            ib = -1;
            ss = 0.;
            ts = 0.;
            iss >> ia >> ib >> ss >> ts;
            ia--;
            ib--;
            pa = numeral(ia,size);
            pb = numeral(ib,size);
            if (ia <= ib) {
               kpi4[npi4] = pa + pb;
            }
            else {
               kpi4[npi4] = pb + pa;
            }
            sslope4[npi4] = ss;
            tslope4[npi4] = ts;
            npi4++;
        }

        // metal ligand field splitting parameters
        else if (keyword == "METAL") {
            ia = -1;
            iss >> ia;
            ia--;
        }

        // biopolymer atom type conversion definitions
        else if (keyword == "BIOTYPE") {
            ia = -1;
            ib = -1;
            iss >> ia;
            ia--;
            getword(record,string,next);
            getstring(record,string,next);
            string = record.substr(next);
            iss.clear();
            iss.str(string);
            iss >> ib;
            ib--;
            if (ia >= maxbio) {
                printf("\n READPRM  --  Too many Biopolymer Types; Increase MAXBIO\n");
                fatal();
            }
            if (ia > -1)  biotyp[ia] = ib;
        }

        // MMFF atom class equivalency parameters
        else if (keyword == "MMFFEQUIV") {
            ia = 999;
            ib = 999;
            ic = 999;
            id = 999;
            ie = 999;
            ig = -1;
            iss >> ia >> ib >> ic >> id >> ie >> ig;
            ia--;
            ib--;
            ic--;
            id--;
            ie--;
            ig--;
            eqclass[0][ig] = ia;
            eqclass[1][ig] = ib;
            eqclass[2][ig] = ic;
            eqclass[3][ig] = id;
            eqclass[4][ig] = ie;
        }

        // MMFF covalent radius and electronegativity parameters
        else if (keyword == "MMFFCOVRAD") {
            ia = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> fc >> bd;
            ia--;
            rad0[ia] = fc;
            paulel[ia] = bd;
        }

        // MMFF atom class property parameters
        else if (keyword == "MMFFPROP") {
            ia = 999;
            ib = 999;
            ic = 999;
            id = 999;
            ie = 999;
            ig = 999;
            ih = 999;
            ii = 999;
            ij = 999;
            iss >> ia >> ib >> ic >> id >> ie >> ig >> ih >> ii >> ij;
            ia--;
            ib--;
            ic--;
            id--;
            ie--;
            ig--;
            ih--;
            ii--;
            ij--;
            crd[ia] = ic;
            val[ia] = id;
            pilp[ia] = ie;
            mltb[ia] = ig;
            arom[ia] = ih;
            lin[ia] = ii;
            sbmb[ia] = ij;
        }

        // MMFF bond stretching parameters
        else if (keyword == "MMFFBOND") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            bt = 2;
            iss >> ia >> ib >> fc >> bd >> bt;
            ia--;
            ib--;
            if (bt == 0) {
               mmff_kb[ib][ia] = fc;
               mmff_kb[ia][ib] = fc;
               mmff_b0[ib][ia] = bd;
               mmff_b0[ia][ib] = bd;
            }
            else if (bt == 1) {
               mmff_kb1[ib][ia] = fc;
               mmff_kb1[ia][ib] = fc;
               mmff_b1[ib][ia] = bd;
               mmff_b1[ia][ib] = bd;
            }
            nb++;
        }

        // MMFF bond stretching empirical rule parameters
        else if (keyword == "MMFFBONDER") {
            ia = -1;
            ib = -1;
            fc = 0.;
            bd = 0.;
            iss >> ia >> ib >> fc >> bd;
            ia--;
            ib--;
            r0ref[ib][ia] = fc;
            r0ref[ia][ib] = fc;
            kbref[ib][ia] = bd;
            kbref[ia][ib] = bd;
        }

        // MMFF bond angle bending parameters
//         else if (keyword == "MMFFANGLE") then
//             ia = 0
//             ib = 0
//             ic = 0
//             fc = 0.0d0
//             an1 = 0.0d0
//             at = 3
//             string = record(next:240)
//             read (string,*,err=680,end=680)  ia,ib,ic,fc,an1,at
//   680       continue
//             na = na + 1
//             if (an1 .ne. 0.0d0) then
//                if (at .eq. 0) then
//                   mmff_ka(ia,ib,ic) = fc
//                   mmff_ka(ic,ib,ia) = fc
//                   mmff_ang0(ia,ib,ic) = an1
//                   mmff_ang0(ic,ib,ia) = an1
//                else if (at .eq. 1) then
//                   mmff_ka1(ia,ib,ic) = fc
//                   mmff_ka1(ic,ib,ia) = fc
//                   mmff_ang1(ia,ib,ic) = an1
//                   mmff_ang1(ic,ib,ia) = an1
//                else if (at .eq. 2) then
//                   mmff_ka2(ia,ib,ic) = fc
//                   mmff_ka2(ic,ib,ia) = fc
//                   mmff_ang2(ia,ib,ic) = an1
//                   mmff_ang2(ic,ib,ia) = an1
//                else if (at .eq. 3) then
//                   mmff_ka3(ia,ib,ic) = fc
//                   mmff_ka3(ic,ib,ia) = fc
//                   mmff_ang3(ia,ib,ic) = an1
//                   mmff_ang3(ic,ib,ia) = an1
//                else if (at .eq. 4) then
//                   mmff_ka4(ia,ib,ic) = fc
//                   mmff_ka4(ic,ib,ia) = fc
//                   mmff_ang4(ia,ib,ic) = an1
//                   mmff_ang4(ic,ib,ia) = an1
//                else if (at .eq. 5) then
//                   mmff_ka5(ia,ib,ic) = fc
//                   mmff_ka5(ic,ib,ia) = fc
//                   mmff_ang5(ia,ib,ic) = an1
//                   mmff_ang5(ic,ib,ia) = an1
//                else if (at .eq. 6) then
//                   mmff_ka6(ia,ib,ic) = fc
//                   mmff_ka6(ic,ib,ia) = fc
//                   mmff_ang6(ia,ib,ic) = an1
//                   mmff_ang6(ic,ib,ia) = an1
//                else if (at .eq. 7) then
//                   mmff_ka7(ia,ib,ic) = fc
//                   mmff_ka7(ic,ib,ia) = fc
//                   mmff_ang7(ia,ib,ic) = an1
//                   mmff_ang7(ic,ib,ia) = an1
//                else if (at .eq. 8) then
//                   mmff_ka8(ia,ib,ic) = fc
//                   mmff_ka8(ic,ib,ia) = fc
//                   mmff_ang8(ia,ib,ic) = an1
//                   mmff_ang8(ic,ib,ia) = an1
//                end if
//             end if
// c
// c     MMFF stretch-bend parameters
// c
//         else if (keyword == "MMFFSTRBND") then
//             ia = 0
//             ib = 0
//             ic = 0
//             abc = 0.0d0
//             cba = 0.0d0
//             sbt = 4
//             string = record(next:240)
//             read (string,*,err=690,end=690)  ia,ib,ic,abc,cba,sbt
//   690       continue
//             if (ia .ne. 0) then
//                if (sbt .eq. 0) then
//                   stbn_abc(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc(ic,ib,ia) = cba
//                   stbn_cba(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba(ic,ib,ia) = abc
//                else if (sbt .eq. 1) then
//                   stbn_abc1(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc1(ic,ib,ia) = cba
//                   stbn_cba1(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba1(ic,ib,ia) = abc
//                else if (sbt .eq. 2) then
//                   stbn_abc2(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc2(ic,ib,ia) = cba
//                   stbn_cba2(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba2(ic,ib,ia) = abc
//                else if (sbt .eq. 3) then
//                   stbn_abc3(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc3(ic,ib,ia) = cba
//                   stbn_cba3(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba3(ic,ib,ia) = abc
//                else if (sbt .eq. 4) then
//                   stbn_abc4(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc4(ic,ib,ia) = cba
//                   stbn_cba4(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba4(ic,ib,ia) = abc
//                else if (sbt .eq. 5) then
//                   stbn_abc5(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc5(ic,ib,ia) = cba
//                   stbn_cba5(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba5(ic,ib,ia) = abc
//                else if (sbt .eq. 6) then
//                   stbn_abc6(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc6(ic,ib,ia) = cba
//                   stbn_cba6(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba6(ic,ib,ia) = abc
//                else if (sbt .eq. 7) then
//                   stbn_abc7(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc7(ic,ib,ia) = cba
//                   stbn_cba7(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba7(ic,ib,ia) = abc
//                else if (sbt .eq. 8) then
//                   stbn_abc8(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc8(ic,ib,ia) = cba
//                   stbn_cba8(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba8(ic,ib,ia) = abc
//                else if (sbt .eq. 9) then
//                   stbn_abc9(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc9(ic,ib,ia) = cba
//                   stbn_cba9(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba9(ic,ib,ia) = abc
//                else if (sbt .eq. 10) then
//                   stbn_abc10(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc10(ic,ib,ia) = cba
//                   stbn_cba10(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba10(ic,ib,ia) = abc
//                else if (sbt .eq. 11) then
//                   stbn_abc11(ia,ib,ic) = abc
//                   if (ic .ne. ia)  stbn_abc11(ic,ib,ia) = cba
//                   stbn_cba11(ia,ib,ic) = cba
//                   if (ic .ne. ia)  stbn_cba11(ic,ib,ia) = abc
//                end if
//             end if
// c
// c     MMFF default stretch-bend parameters
// c
//         else if (keyword == "MMFFDEFSTBN") then
//             string = record(next:240)
//             ia = 1000
//             ib = 1000
//             ic = 1000
//             abc = 0.0d0
//             cba = 0.0d0
//             read (string,*,err=700,end=700)  ia,ib,ic,abc,cba
//   700       continue
//             defstbn_abc(ia,ib,ic) = abc
//             defstbn_cba(ia,ib,ic) = cba
//             defstbn_abc(ic,ib,ia) = cba
//             defstbn_cba(ic,ib,ia) = abc
// c
// c     MMFF out-of-plane bend parameters
// c
//         else if (keyword == "MMFFOPBEND") then
//             ia = 0
//             ib = 0
//             ic = 0
//             id = 0
//             fc = 0.0d0
//             string = record(next:240)
//             read (string,*,err=710,end=710)  ia,ib,ic,id,fc
//   710       continue
//             call numeral (ia,pa,size)
//             call numeral (ib,pb,size)
//             call numeral (ic,pc,size)
//             call numeral (id,pd,size)
//             nopb = nopb + 1
//             if (ic .le. id) then
//                kopb(nopb) = pa//pb//pc//pd
//             else
//                kopb(nopb) = pa//pb//pd//pc
//             end if
//             opbn(nopb) = fc
// c           if (ic.gt.0 .or. id.gt.0) then
// c              nopb = nopb + 1
// c              if (ib .le. id) then
// c                 kopb(nopb) = pc//pb//pb//pd
// c              else
// c                 kopb(nopb) = pc//pb//pd//pb
// c              end if
// c              opbn(nopb) = fc
// c              nopb = nopb + 1
// c              if (ia .le. ic) then
// c                 kopb(nopb) = pd//pb//pa//pc
// c              else
// c                 kopb(nopb) = pd//pb//pc//pa
// c              end if
// c              opbn(nopb) = fc
// c           end if
// c
// c     MMFF torsional parameters
// c
//         else if (keyword == "MMFFTORSION") then
//             ia = 0
//             ib = 0
//             ic = 0
//             id = 0
//             do i = 1, 6
//                vt(i) = 0.0d0
//                st(i) = 0.0d0
//                ft(i) = 0
//             end do
//             tt = 3
//             string = record(next:240)
//             read (string,*,err=720,end=720)  ia,ib,ic,id,(vt(j),
//      &                                       st(j),ft(j),j=1,3),tt
//   720       continue
//             call numeral (ia,pa,size)
//             call numeral (ib,pb,size)
//             call numeral (ic,pc,size)
//             call numeral (id,pd,size)
//             nt = nt + 1
//             if (tt .eq. 0) then
//                if (ib .lt. ic) then
//                   kt(nt) = pa//pb//pc//pd
//                else if (ic .lt. ib) then
//                   kt(nt) = pd//pc//pb//pa
//                else if (ia .le. id) then
//                   kt(nt) = pa//pb//pc//pd
//                else if (id .lt. ia) then
//                   kt(nt) = pd//pc//pb//pa
//                end if
//                call torphase (ft,vt,st)
//                t1(1,nt) = vt(1)
//                t1(2,nt) = st(1)
//                t2(1,nt) = vt(2)
//                t2(2,nt) = st(2)
//                t3(1,nt) = vt(3)
//                t3(2,nt) = st(3)
//             else if (tt .eq. 1) then
//                if (ib .lt. ic) then
//                   kt_1(nt) = pa//pb//pc//pd
//                else if (ic .lt. ib) then
//                   kt_1(nt) = pd//pc//pb//pa
//                else if (ia .le. id) then
//                   kt_1(nt) = pa//pb//pc//pd
//                else if (id .lt. ia) then
//                   kt_1(nt) = pd//pc//pb//pa
//                end if
//                call torphase (ft,vt,st)
//                t1_1(1,nt) = vt(1)
//                t1_1(2,nt) = st(1)
//                t2_1(1,nt) = vt(2)
//                t2_1(2,nt) = st(2)
//                t3_1(1,nt) = vt(3)
//                t3_1(2,nt) = st(3)
//             else if (tt .eq. 2) then
//                if (ib .lt. ic) then
//                   kt_2(nt) = pa//pb//pc//pd
//                else if (ic .lt. ib) then
//                   kt_2(nt) = pd//pc//pb//pa
//                else if (ia .le. id) then
//                   kt_2(nt) = pa//pb//pc//pd
//                else if (id .lt. ia) then
//                   kt_2(nt) = pd//pc//pb//pa
//                end if
//                call torphase (ft,vt,st)
//                t1_2(1,nt) = vt(1)
//                t1_2(2,nt) = st(1)
//                t2_2(1,nt) = vt(2)
//                t2_2(2,nt) = st(2)
//                t3_2(1,nt) = vt(3)
//                t3_2(2,nt) = st(3)
//             else if (tt .eq. 4) then
//                nt4 = nt4 + 1
//                if (ib .lt. ic) then
//                   kt4(nt4) = pa//pb//pc//pd
//                else if (ic .lt. ib) then
//                   kt4(nt4) = pd//pc//pb//pa
//                else if (ia .le. id) then
//                   kt4(nt4) = pa//pb//pc//pd
//                else if (id .lt. ia) then
//                   kt4(nt4) = pd//pc//pb//pa
//                end if
//                call torphase (ft,vt,st)
//                t14(1,nt4) = vt(1)
//                t14(2,nt4) = st(1)
//                t24(1,nt4) = vt(2)
//                t24(2,nt4) = st(2)
//                t34(1,nt4) = vt(3)
//                t34(2,nt4) = st(3)
//             else if (tt .eq. 5) then
//                nt5 = nt5 + 1
//                if (ib .lt. ic) then
//                   kt5(nt5) = pa//pb//pc//pd
//                else if (ic .lt. ib) then
//                   kt5(nt5) = pd//pc//pb//pa
//                else if (ia .le. id) then
//                   kt5(nt5) = pa//pb//pc//pd
//                else if (id .lt. ia) then
//                   kt5(nt5) = pd//pc//pb//pa
//                end if
//                call torphase (ft,vt,st)
//                t15(1,nt5) = vt(1)
//                t15(2,nt5) = st(1)
//                t25(1,nt5) = vt(2)
//                t25(2,nt5) = st(2)
//                t35(1,nt5) = vt(3)
//                t35(2,nt5) = st(3)
//             end if
// c
// c     MMFF van der Waals parameters
// c
//         else if (keyword == "MMFFVDW") then
//             ia = 0
//             rd = 0.0d0
//             ep = 0.0d0
//             rdn = 0.0d0
//             da1 = 'C'
//             string = record(next:240)
//             read (string,*,err=730,end=730)  ia,rd,alphi,nni,gi,da1
//   730       continue
//             if (ia .ne. 0) then
//                rad(ia) = rd
//                g(ia) = gi
//                alph(ia) = alphi
//                nn(ia) = nni
//                da(ia) = da1
//             end if
// c
// c     MMFF bond charge increment parameters
// c
//         else if (keyword == "MMFFBCI") then
//             ia = 0
//             ib = 0
//             cg = 1000.0d0
//             bt = 2
//             string = record(next:240)
//             read (string,*,err=740,end=740)  ia,ib,cg,bt
//   740       continue
//             if (ia .ne. 0) then
//                if (bt .eq. 0) then
//                   bci(ia,ib) = cg
//                   bci(ib,ia) = -cg
//                else if (bt .eq. 1) then
//                   bci_1(ia,ib) = cg
//                   bci_1(ib,ia) = -cg
//                end if
//             end if
// c
// c     MMFF partial bond charge increment parameters
// c
//         else if (keyword == "MMFFPBCI") then
//             ia = 0
//             string = record(next:240)
//             read (string,*,err=750,end=750)  ia,cg,factor
//   750       continue
//             if (ia .ne. 0) then
//                pbci(ia) = cg
//                fcadj(ia) = factor
//             end if
// c
// c     MMFF aromatic ion parameters
// c
//         else if (keyword == "MMFFAROM") then
//             string = record(next:240)
//             read (string,*,err=760,end=760)  ia,ib,ic,id,ie,if
//   760       continue
//             if (ie.eq.0 .and. id.eq.0) then
//                mmffarom(ia,if) = ic
//             else if (id .eq. 1) then
//                mmffaromc(ia,if) = ic
//             else if (ie .eq. 1) then
//                mmffaroma(ia,if) = ic
//             end if
//          end if
//       end do
        iprm++;
    }
}
