/* ====================================================================

This file contains a series of procedures used to compute the
intrinsic volumes of a union of balls, i.e. its surface area, volume,
mean integrated curvature, and Gaussian curvature, and optionally their derivatives
with respect to the coordinates of the centers of the ball.

 ==================================================================== */

#pragma once

#include "edge.h"
#include "face.h"
#include "gausscorners.h"
#include "tetra.h"
#include "tetrahedron.h"
#include "vector.h"
#include <vector>
 
/* ====================================================================
   class
 ==================================================================== */

class VOLUMES {
public:

    void ball_dvolumes(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
        std::vector<Edge>& edges, std::vector<Face>& faces,
        double *WSurf, double *WVol, double *WMean, double *WGauss,
        double *Surf, double *Vol, double *Mean, double *Gauss,
        double *ballwsurf, double *ballwvol, double *ballwmean, double *ballwgauss,
        double *dsurf_coord, double *dvol_coord, double *dmean_coord,
        double *dgauss_coord, int option);

private:

    double distance2(std::vector<Vertex>& vertices, int n1, int n2);

    void twosphere_info(double ra, double ra2, double rb, double rb2,
        double rab, double rab2, double *surfa, double *surfb,
        double *vola, double *volb, double *r, double *phi, double *l);

    void twosphere_dinfo(double ra, double ra2, double rb, double rb2,
        double rab, double rab2, double *surfa, double *surfb,
        double *vola, double *volb, double *r, double *phi, double *l,
        double *dsurfa, double *dsurfb, double *dvola, double *dvolb, 
        double *dr, double *dphi, double *dl, int option);

    void threesphere_dvol(double ra, double rb,double rc, double ra2,
        double rb2, double rc2, double rab, double rac, double rbc,
        double rab2, double rac2, double rbc2, double *angle, double deriv[6][3],
        double *surfa, double *surfb, double *surfc, double *vola, double *valb,
        double *volc, double *dsurfa, double *dsurfb, double *dsurfc, 
        double *dvola, double *dvolb, double *dvolc, int option);

    double plane_dist(double ra2, double rb2, double rab2);

    double pi = M_PI;
    double twopi = 2.0*pi;

    double coef_E = 1.;
    double coef_F = 2.;
};
