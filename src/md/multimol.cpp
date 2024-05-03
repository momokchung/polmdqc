// Author: Moses KJ Chung
// Year:   2024

#include "alfboxsize.h"
#include "alphamol.h"
#include "alphmol.h"
#include "dlauny2.h"
#include "hilbert.h"
#include "multimol.h"
#include <pthread.h>

namespace polmdqc
{
#define NUM_THREADS 128 
pthread_t threads[NUM_THREADS];

typedef struct thread_data {
    bool deriv;
    int N1;
    int N2;
    int nlist2;
    real buffer;
    real wsurf;
    real wvol;
    real wmean;
    real wgauss;
    real* surf;
    real* vol;
    real* mean;
    real* gauss;
    real* dsurf;
    real* dvol;
    real* dmean;
    real* dgauss;
} thread_data;

thread_data info[NUM_THREADS];
int threadids[NUM_THREADS];

inline bool inBox(real Point[3], real Xlim[2], real Ylim[2], real Zlim[2])
{
    if(Point[0] < Xlim[0] || Point[0] >= Xlim[1]) return false;
    if(Point[1] < Ylim[0] || Point[1] >= Ylim[1]) return false;
    if(Point[2] < Zlim[0] || Point[2] >= Zlim[1]) return false;
    return true;
}

/////////////////////////////////////////////////////
//                                                 //
//  singlemol  --  launch singlethreaded AlphaMol  //
//                                                 //
/////////////////////////////////////////////////////

void* singlemol(void* data)
{
    real xmin,ymin,zmin;
    real xmax,ymax,zmax;
    real r;
    real Xbox_buf[2],Ybox_buf[2],Zbox_buf[2];

    int threadid = *((int *) data);
    int N1       = info[threadid].N1;
    int N2       = info[threadid].N2;
    int natm     = N2 - N1;
    int natoms   = alfatoms.size();
    real buffer  = info[threadid].buffer;

    alfboxsize(&alfatoms[N1], natm, xmin, ymin, zmin, xmax, ymax, zmax, r);

    Xbox_buf[0] = xmin - buffer; Xbox_buf[1] = xmax + buffer;
    Ybox_buf[0] = ymin - buffer; Ybox_buf[1] = ymax + buffer;
    Zbox_buf[0] = zmin - buffer; Zbox_buf[1] = zmax + buffer;

    int nlist2;
    int *list2 = new int[natoms];
    real Point[3];
    bool test;

    nlist2 = 0;
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < 3; j++) Point[j] = alfatoms[i].coord[j];
        test = inBox(Point, Xbox_buf, Ybox_buf, Zbox_buf);
        if (test) {
            list2[nlist2] = i;
            nlist2++;
        }
    }
    for (int i = N2; i < natoms; i++) {
        for (int j = 0; j < 3; j++) Point[j] = alfatoms[i].coord[j];
        test = inBox(Point, Xbox_buf, Ybox_buf, Zbox_buf);
        if (test) {
            list2[nlist2] = i;
            nlist2++;
        }
    }

    int ntot = natm + nlist2;

    AlfAtom *newatoms = new AlfAtom[ntot];

    for (int i = 0; i < natm; i++) {
        newatoms[i] = alfatoms[N1+i];
        // newatoms[i].index = i;
    }
    // TODO: might not need
    sort3DHilbert(newatoms, natm, 0, 0, xmax, ymax, zmax, xmax, ymax, zmax, 0);

    int nat = natm;
    for (int i = 0; i < nlist2; i++) {
        int k = list2[i];
        newatoms[nat] = alfatoms[k];
        nat++;
    }

    int nfudge = 8;
    real tmp1;
    real tmp2;
    real tmp3;
    real tmp4;
    real* surfthd = new real[ntot+nfudge];
    real* volthd = new real[ntot+nfudge];
    real* meanthd = new real[ntot+nfudge];
    real* gaussthd = new real[ntot+nfudge];
    real* dsurfthd = new real[3*(ntot+nfudge)];
    real* dvolthd = new real[3*(ntot+nfudge)];
    real* dmeanthd = new real[3*(ntot+nfudge)];
    real* dgaussthd = new real[3*(ntot+nfudge)];
    bool deriv = info[threadid].deriv;

    alphamol(ntot, newatoms, tmp1, tmp2, tmp3, tmp4,
        surfthd, volthd, meanthd, gaussthd,
        dsurfthd,dvolthd, dmeanthd, dgaussthd, deriv);

    // transfer information to thread
    info[threadid].wsurf = 0;
    info[threadid].wvol = 0;
    info[threadid].wmean = 0;
    info[threadid].wgauss = 0;

    for(int i = 0; i < natm; i++) {
        info[threadid].wsurf += surfthd[i];
        info[threadid].wvol += volthd[i];
        info[threadid].wmean += meanthd[i];
        info[threadid].wgauss += gaussthd[i];
    }

    for (int i = 0; i < natm; i++) {
        info[threadid].surf[newatoms[i].index] = surfthd[i];
        info[threadid].vol[newatoms[i].index] = volthd[i];
        info[threadid].mean[newatoms[i].index] = meanthd[i];
        info[threadid].gauss[newatoms[i].index] = gaussthd[i];
        if (deriv) {
            for (int j = 0; j < 3; j++) {
                info[threadid].dsurf[3*newatoms[i].index + j] = dsurfthd[3*i+j];
                info[threadid].dvol[3*newatoms[i].index + j] = dvolthd[3*i+j];
                info[threadid].dmean[3*newatoms[i].index + j] = dmeanthd[3*i+j];
                info[threadid].dgauss[3*newatoms[i].index + j] = dgaussthd[3*i+j];
            }
        }
    }

    // info[threadid].nlist2 = nlist2;

    delete[] newatoms;
    delete[] surfthd;
    delete[] volthd;
    delete[] meanthd;
    delete[] gaussthd;
    delete[] dsurfthd;
    delete[] dvolthd;
    delete[] dmeanthd;
    delete[] dgaussthd;
    delete[] list2;

    return 0;
}

///////////////////////////////////////////////////
//                                               //
//  multimol  --  launch multithreaded AlphaMol  //
//                                               //
///////////////////////////////////////////////////

// "multimol" launches a multithreaded AlphaMol

void multimol(real buffer, bool deriv, int nthreads, std::vector<int>& Nval)
{
    int N1,N2;
    int natoms = alfatoms.size();
    int nval = natoms/nthreads;

    for (int i = 0; i < nthreads; i++) {
        if (Nval[nthreads] == 0) {
            N1 = i*nval;
            N2 = N1 + nval;
            if (i == nthreads-1) N2 = natoms;
        }
        else {
            N1 = Nval[i];
            N2 = Nval[i+1];
        }

        threadids[i] = i;

        info[i].deriv  = deriv;
        info[i].N1     = N1;
        info[i].N2     = N2;
        info[i].buffer = buffer;
        info[i].surf   = surf.ptr();
        info[i].vol    = vol.ptr();
        info[i].mean   = mean.ptr();
        info[i].gauss  = gauss.ptr();
        info[i].dsurf  = dsurf.ptr();
        info[i].dvol   = dvol.ptr();
        info[i].dmean  = dmean.ptr();
        info[i].dgauss = dgauss.ptr();

        pthread_create(&threads[i], NULL, singlemol, (void*) &threadids[i]);
    }

    // for(int i = 0; i < natoms; i++) {
    //     std::cout << " i: " << i << " atom: " << alfatoms[i].index << std::endl;
    // }

    wsurf = 0;
    wvol = 0;
    wmean = 0;
    wgauss = 0;
    for (int i = 0; i < nthreads; i++) {
        pthread_join(threads[i], NULL);
        wsurf += info[i].wsurf;
        wvol += info[i].wvol;
        wmean += info[i].wmean;
        wgauss += info[i].wgauss;
    }
}
}
