// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alphamol.h"
#include "alphmol.h"
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
// 
// literature reference:
// 
// P. Koehl, A. Akopyan, and H. Edelsbrunner, "Computing the Volume,
// Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
// Derivatives", Journal of Chemical Information and Modeling,
// 63, 973-985, (2023).
// 
// github reference:
// 
// https://github.com/pkoehl/AlphaMol

void alphamol(double r_h2o, bool computeDeriv)
{
    DELCX delcx;
    ALFCX alfcx;
    VOLUMES volumes;

    int flag_deriv = 0;
    if (computeDeriv) flag_deriv = 1;

    // // print information about system
    // std::cout << "" << std::endl;
	// std::cout << " Input file                : " << filename << std::endl;
	// std::cout << " Number of atoms (balls)   : " << n << std::endl;
	// std::cout << " Probe radius              : " << r_h2o << std::endl;
	// std::cout << "" << std::endl;

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
	// std::cout << "Delaunay compute time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

    // Generate alpha complex (with alpha=0.0)
	start_s = clock();
	double alpha = 0;
	alfcx.alfcx(alpha, vertices, tetra);
	stop_s = clock();
	// std::cout << "AlphaCx compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
	std::vector<Edge> edges;
	std::vector<Face> faces;
	alfcx.alphacxEdges(tetra, edges);
	alfcx.alphacxFaces(tetra, faces);

	int nfudge = 8;
    surf.allocate(natoms+nfudge);
    if (computeDeriv) {
        dsurf.allocate(3*(natoms+nfudge));
        memset(dsurf.ptr(), 0, 3*(natoms+nfudge)*sizeof(double));
    }

	vol.allocate(natoms+nfudge);
    if (computeDeriv) {
        dvol.allocate(3*(natoms+nfudge));
        memset(dvol.ptr(), 0, 3*(natoms+nfudge)*sizeof(double));
    }

	mean.allocate(natoms+nfudge);
    if (computeDeriv) {
        dmean.allocate(3*(natoms+nfudge));
        memset(dmean.ptr(), 0, 3*(natoms+nfudge)*sizeof(double));
    }

	gauss.allocate(natoms+nfudge);
    if (computeDeriv) {
        dgauss.allocate(3*(natoms+nfudge));
        memset(dgauss.ptr(), 0, 3*(natoms+nfudge)*sizeof(double));
    }

	start_s = clock();
	volumes.ball_dvolumes(vertices, tetra, edges, faces, &wsurf, &wvol, &wmean, &wgauss, &tsurf, &tvol, &tmean, &tgauss,
	    surf.ptr(), vol.ptr(), mean.ptr(), gauss.ptr(), dsurf.ptr(), dvol.ptr(), dmean.ptr(), dgauss.ptr(), flag_deriv);
	stop_s = clock();

	// std::cout << "Volumes compute time : " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds" << std::endl;

	// std::cout << " " << std::endl;
	// std::cout << "Biomolecule from file      : " << filename << std::endl;
	// std::cout << "Number of atoms (balls)    : " << natoms << std::endl;
	// std::cout << "Probe radius               : " << r_h2o << std::endl;
	// std::cout << "Unweighted surface area    : " << std::setw(16) << std::fixed << std::setprecision(8) << tsurf << std::endl;
	// std::cout << "Weighted surface area      : " << std::setw(16) << std::fixed << std::setprecision(8) << wsurf << std::endl;
	// std::cout << "Unweighted volume          : " << std::setw(16) << std::fixed << std::setprecision(8) << tvol << std::endl;
	// std::cout << "Weighted volume            : " << std::setw(16) << std::fixed << std::setprecision(8) << wvol << std::endl;
	// std::cout << "Unweighted mean curvature  : " << std::setw(16) << std::fixed << std::setprecision(8) << tmean << std::endl;
	// std::cout << "Weighted mean curvature    : " << std::setw(16) << std::fixed << std::setprecision(8) << wmean << std::endl;
	// std::cout << "Unweighted Gauss curvature : " << std::setw(16) << std::fixed << std::setprecision(8) << tgauss << std::endl;
	// std::cout << "Weighted Gauss curvature   : " << std::setw(16) << std::fixed << std::setprecision(8) << wgauss << std::endl;
	// std::cout << " " << std::endl;

    // std::cout << "Surface derivatives: " << std::endl;
	// for(int i = 0; i < natoms; i++) {
    //     printf("%24.16e%24.16e%24.16e\n", dsurf[3*i], dsurf[3*i+1], dsurf[3*i+2]);
	// }
	// std::cout << " " << std::endl;

    // std::cout << "Volume derivatives: " << std::endl;
	// for(int i = 0; i < natoms; i++) {
    //     printf("%24.16e%24.16e%24.16e\n", dvol[3*i], dvol[3*i+1], dvol[3*i+2]);
	// }
	// std::cout << " " << std::endl;

	// std::cout << "Mean curvature derivatives: " << std::endl;
	// for(int i = 0; i < natoms; i++) {
    //     printf("%24.16e%24.16e%24.16e\n", dmean[3*i], dmean[3*i+1], dmean[3*i+2]);
	// }
	// std::cout << " " << std::endl;

	// std::cout << "Gauss curvature derivatives: " << std::endl;
	// for(int i = 0; i < natoms; i++) {
    //     printf("%24.16e%24.16e%24.16e\n", dgauss[3*i], dgauss[3*i+1], dgauss[3*i+2]);
	// }
	// std::cout << " " << std::endl;

	delete [] coord;
    delete [] radii;
    delete [] coefS;
    delete [] coefV;
    delete [] coefM;
    delete [] coefG;
}
}
