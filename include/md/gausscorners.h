/* ====================================================================
	sets of tools to compute the contributions of the corners
	of the space-filling diagram to the Gaussian curvature
 ==================================================================== */

#pragma once

#include <cmath>

/* ====================================================================
   class
 ==================================================================== */

class GAUSSCORNER {

public:

    void threesphere_dgauss(double ra, double rb, double rc, 
    double ra2, double rb2, double rc2,
    double rab, double rac, double rbc,
    double rab2, double rac2, double rbc2,
    double *areaA, double *areaB, double *areaC, double darea[3][3], int option);

private:

    double trig_darea(double a, double b, double c, double *der_S, int option);

    double trig_dradius(double a, double b, double c, double *der_r, int option);

    double sign(double a, double b, double c);
};
