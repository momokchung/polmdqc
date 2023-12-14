/* ==============================================================================
 *	Sos_gmp
 *									  
 *  Performs all operations   
 *  with multi precision arithmetics, using the package GMP               
 *									  
   ============================================================================== */

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"

#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

/* ==============================================================================
   The SOS class
   ============================================================================== */

class SOS 
{
public:

	void real_to_gmp(double coord, mpz_t val);

	void build_weight(mpz_t ax, mpz_t ay, mpz_t az, mpz_t r, mpz_t w);

	void sos_minor2_gmp(double xa, double xb, int *res);

	void sos_minor3_gmp(double xa, double ya, double xb, double yb, double xc, double yc, int *res);

	void sos_minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res);

	void sos_minor5_gmp(double *coord_a, double ra, double *coord_b, double rb, 
	double *coord_c, double rc, double *coord_d, double rd, double *coord_e, double re,
	int *res) ;

	void minor4_gmp(double *coord_a, double *coord_b, double *coord_c, double *coord_d, int *res);

	void init_sos_gmp();

	void clear_sos_gmp();

private:

	void deter2_gmp(mpz_t deter, mpz_t a, mpz_t b);

	void deter3_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a21, mpz_t a22,
		mpz_t a31, mpz_t a32);

	void deter4_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a21,
		mpz_t a22, mpz_t a23, mpz_t a31, mpz_t a32,
		mpz_t a33, mpz_t a41, mpz_t a42, mpz_t a43);

	void deter5_gmp(mpz_t deter, mpz_t a11, mpz_t a12, mpz_t a13, mpz_t a14,
		mpz_t a21, mpz_t a22, mpz_t a23, mpz_t a24,
		mpz_t a31, mpz_t a32, mpz_t a33, mpz_t a34,
		mpz_t a41, mpz_t a42, mpz_t a43, mpz_t a44,
		mpz_t a51, mpz_t a52, mpz_t a53, mpz_t a54);

	mpz_t a11_mp,a12_mp,a13_mp,a14_mp;
	mpz_t a21_mp,a22_mp,a23_mp,a24_mp;
	mpz_t a31_mp,a32_mp,a33_mp,a34_mp;
	mpz_t a41_mp,a42_mp,a43_mp,a44_mp;
	mpz_t a51_mp,a52_mp,a53_mp,a54_mp;
	mpz_t r1_mp,r2_mp, r3_mp, r4_mp, r5_mp;

	mpz_t temp1,temp2,temp3,temp4;
	mpz_t val1,val2,val3;

	mpz_t c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,c41,c42,c43,c44;
	mpz_t d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3;

	double scale;
};
