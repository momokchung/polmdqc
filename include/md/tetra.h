/* ====================================================================
	sets of tools to measure tetrahedron
 ==================================================================== */

#pragma once

#include <cmath>

/* ====================================================================
   class
 ==================================================================== */

class TETRAGEOM {
public:

    void tetra_dihed(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq, double r34sq, double *angle,
    double *cosine, double *sine);

    void tetra_dihed_der(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq, double r34sq, double *angle,
    double *cosine, double *sine, double deriv[6][6]);

    void tetra_dihed_der3(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq, double r34sq, double *angle,
    double *cosine, double *sine, double deriv[6][3], int option);

    void tetra_3dihed_cos(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq,double r34sq, double *cosine);

    void tetra_3dihed_dcos(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq,double r34sq, double *cosine, 
    double deriv[3][3], int option);

    double tetra_volume(double r12sq, double r13sq, double r14sq,
    double r23sq, double r24sq, double r34sq);

    void tetra_Voronoi(double ra2,double rb2,double rc2,double rd2,
        double rab, double rac, double rad, double rbc, double rbd,
        double rcd, double rab2, double rac2, double rad2,double rbc2,
        double rbd2, double rcd2, double *cos_ang, double *sin_ang,
        double *vola, double *volb, double *volc, double *vold);

    void tetra_Voronoi_der(double ra2,double rb2,double rc2,double rd2,
        double rab, double rac, double rad, double rbc, double rbd,
        double rcd, double rab2, double rac2, double rad2,double rbc2,
        double rbd2, double rcd2, double *cos_ang, double *sin_ang,
        double deriv[6][6], double *vola, double *volb, double *volc, 
        double *vold, double *dvola, double *dvolb, double *dvolc,
        double *dvold, int option);

    double plane_dist(double ra2, double rb2, double rab2);

private:

    double pi = M_PI;
    double twopi = 2.0*pi;
};
