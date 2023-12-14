// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alphamol.h"
#include "atomid.h"
#include "atoms.h"
#include "files.h"
#include "kvdws.h"
#include <iostream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

// "alphamol" computes volume, surface area, mean, and gaussian curvature

// literature reference:

// P. Koehl, A. Akopyan, and H. Edelsbrunner, "Computing the Volume,
// Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
// Derivatives", Journal of Chemical Information and Modeling,
// 63, 973-985, (2023).

// github reference:

// https://github.com/pkoehl/AlphaMol

void alphamol(double r_h2o, int flag_deriv)
{
    DELCX delcx;
    ALFCX alfcx;
    VOLUMES volumes;

    // print information about system
    std::cout << "" << std::endl;
	std::cout << " Input file                : " << filename << std::endl;
	std::cout << " Number of atoms (balls)   : " << n << std::endl;
	std::cout << " Probe radius              : " << r_h2o << std::endl;
	std::cout << "" << std::endl;

    // Compute Delaunay triangulation
	clock_t start_s, stop_s;

	std::vector<Vertex> vertices;
	std::vector<Tetrahedron> tetra;

	int natoms = n;
	double *coord = new double[3*natoms];
	double *radii = new double[natoms];
	double *coefS = new double[natoms];
	double *coefV = new double[natoms];
	double *coefM = new double[natoms];
	double *coefG = new double[natoms];

	for(int i = 0; i < natoms; i++) {
        coord[3*i] = x[i];
        coord[3*i+1] = y[i];
        coord[3*i+2] = z[i];
        radii[i] = rad[atomClass[i]] + r_h2o;

		coefS[i] = 1.0;
		coefV[i] = 1.0;
		coefM[i] = 1.0;
		coefG[i] = 1.0;
	}
	
	delcx.setup(natoms, coord, radii, coefS, coefV, coefM, coefG, vertices, tetra);

	start_s = clock();
	delcx.regular3D(vertices, tetra);
	stop_s = clock();
	std::cout << "Delaunay compute time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

    // Generate alpha complex (with alpha=0.0)
	start_s = clock();
	double alpha = 0;
	alfcx.alfcx(alpha, vertices, tetra);
	stop_s = clock();
	std::cout << "AlphaCx compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
	std::vector<Edge> edges;
	std::vector<Face> faces;
	alfcx.alphacxEdges(tetra, edges);
	alfcx.alphacxFaces(tetra, faces);

	double Surf, WSurf, Vol, WVol, Mean, WMean, Gauss, WGauss;

	int nfudge = 8;
	double *ballwsurf = new double[natoms+nfudge];
	double *dsurf = new double[3*(natoms+nfudge)];
	memset(dsurf, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwvol, *dvol;
	ballwvol = new double[natoms+nfudge];
	dvol = new double[3*(natoms+nfudge)];
	memset(dvol, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwmean, *dmean;
	ballwmean = new double[natoms+nfudge];
	dmean = new double[3*(natoms+nfudge)];
	memset(dmean, 0, 3*(natoms+nfudge)*sizeof(double));

	double *ballwgauss, *dgauss;
	ballwgauss = new double[natoms+nfudge];
	dgauss = new double[3*(natoms+nfudge)];
	memset(dgauss, 0, 3*(natoms+nfudge)*sizeof(double));

	start_s = clock();
	volumes.ball_dvolumes(vertices, tetra, edges, faces, &WSurf, &WVol,
	&WMean, &WGauss, &Surf, &Vol, &Mean, &Gauss, ballwsurf, ballwvol,
	ballwmean, ballwgauss, dsurf, dvol, dmean, dgauss, flag_deriv);
	stop_s = clock();

	std::cout << "Volumes compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

	std::cout << " " << std::endl;
	std::cout << "Biomolecule from file      : " << filename << std::endl;
	std::cout << "Number of atoms (balls)    : " << natoms << std::endl;
	std::cout << "Probe radius               : " << r_h2o << std::endl;
	std::cout << "Unweighted surface area    : " << std::setw(16) << std::fixed << std::setprecision(8) << Surf << std::endl;
	std::cout << "Weighted surface area      : " << std::setw(16) << std::fixed << std::setprecision(8) << WSurf << std::endl;
	std::cout << "Unweighted volume          : " << std::setw(16) << std::fixed << std::setprecision(8) << Vol << std::endl;
	std::cout << "Weighted volume            : " << std::setw(16) << std::fixed << std::setprecision(8) << WVol << std::endl;
	std::cout << "Unweighted mean curvature  : " << std::setw(16) << std::fixed << std::setprecision(8) << Mean << std::endl;
	std::cout << "Weighted mean curvature    : " << std::setw(16) << std::fixed << std::setprecision(8) << WMean << std::endl;
	std::cout << "Unweighted Gauss curvature : " << std::setw(16) << std::fixed << std::setprecision(8) << Gauss << std::endl;
	std::cout << "Weighted Gauss curvature   : " << std::setw(16) << std::fixed << std::setprecision(8) << WGauss << std::endl;
	std::cout << " " << std::endl;

	delete [] coord; delete [] radii; delete [] coefS; delete [] coefV; delete [] coefM; delete [] coefG;
	delete [] ballwsurf; delete [] dsurf; delete [] ballwvol; delete [] dvol;
	delete [] ballwmean; delete [] dmean; delete [] ballwgauss; delete [] dgauss;
}
}
