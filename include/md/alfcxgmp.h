/* =======================================================================================
	alfx_gmp

        Performs all predicates for alfcx using gmp
                                                                        
 ======================================================================================= */

#pragma once

#include "gmp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =======================================================================================
	the class
 ======================================================================================= */

class ALFCX_GMP {

public:

    void tetra_radius_gmp(double *a, double *b, double *c, double *d, double ra,
    double rb, double rc, double rd, int *test, double alpha);

    void vertex_attach_gmp(double *a, double *b, double ra, double rb, int *testa, int *testb);

    void edge_attach_gmp(double *a, double *b, double *c, double ra,
    double rb, double rc, int *test, int *memory);

    void edge_radius_gmp(double *a, double *b, double ra, double rb,
    int *test, double alpha, int *memory);

    void triangle_attach_gmp(double *a, double *b, double *c, double *d,
        double ra, double rb, double rc, double rd, int *test, int *memory);

    void triangle_radius_gmp(double *a, double *b, double *c, double ra,
    double rb, double rc, int *test, double alpha, int *memory);

    void set_alf_gmp();
    void clear_alf_gmp();

private:

    void set_edge(double *a, double *b, double ra, double rb);

    void set_triangle(double *a, double *b, double *c, double ra, double rb, double rc);

    void real_to_gmp(double *coord, int idx, mpz_t val);

    void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

    void scalar_to_gmp(double coord, mpz_t val);

    mpz_t temp1, temp2, temp3;

    mpz_t ra2,rb2,dist2,dtest, num, den;
    mpz_t r_11, r_22, r_33, r_14, r_313, r_212,diff, det0, det1, det2, det3, det4;
    mpz_t Dabc, Dabd, Dacd, Dbcd, Dabcd;
    mpz_t wa,wb,wc,wd;
    
    mpz_t ra_mp,rb_mp, rc_mp, rd_mp;
    mpz_t alp;
    
    mpz_t res[4][5], res2_c[4][5];
    mpz_t a_mp[5], b_mp[5], c_mp[5], d_mp[5];
    mpz_t Tab[4], Sab[4], Dab[5];
    mpz_t Sac[4], Sad[4], Sbc[4], Sbd[4], Scd[4];
    mpz_t Sa[4], Sb[4], Sd[4];
    mpz_t Sam1[4], Sbm1[4], Scm1[4], Sdm1[4];
    mpz_t Deter[4];
    mpz_t Tc[4],Sc[4];
    mpz_t Mab[4][5], Mac[4][5], Mbc[4][5], S[4][5], T[4][5];

    double scale = 1.e8;
};
