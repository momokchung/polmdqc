/* -------------------------------------------------------------------------- *
 *                                 GaussVol                                   *
 * -------------------------------------------------------------------------- *
 * This file is part of the AGBNP/OpenMM implicit solvent model software      *
 * implementation funded by the National Science Foundation under grant:      *
 * NSF SI2 1440665  "SI2-SSE: High-Performance Software for Large-Scale       *
 * Modeling of Binding Equilibria"                                            *
 *                                                                            *
 * copyright (c) 2016 Emilio Gallicchio                                       *
 * Authors: Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>                 *
 * Contributors:                                                              *
 *                                                                            *
 *  AGBNP/OpenMM is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Lesser General Public License version 3     *
 *  as published by the Free Software Foundation.                             *
 *                                                                            *
 *  AGBNP/OpenMM is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with AGBNP/OpenMM.  If not, see <http://www.gnu.org/licenses/>      *
 *                                                                            *
 * -------------------------------------------------------------------------- */

#pragma once

#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// conversion factors from spheres to Gaussians
constexpr double kfc = 2.2269859253;

// have switching function
constexpr double min_gvol = std::numeric_limits<double>::min();

// maximum overlap level
constexpr int max_order = 16;

//use nm and kj
constexpr double ang3 = 1.;

//volume cutoffs in switching function
constexpr double volmina = 0.01 * ang3;
constexpr double volminb = 0.1 * ang3;

// 3D Gaussian, V,c,a representation
class GaussianVca {
public:
    double v;  // Gaussian volume
    double a;  // Gaussian exponent
    double cx; // center x-coordinate
    double cy; // center y-coordinate
    double cz; // center z-coordinate
};

// switching function used in Gaussian overlap function
double pol_switchfunc(double gvol, double volmina, double volminb, double& sp);

// overlap between two Gaussians represented by a (V,c,a) triplet
// V: volume of Gaussian
// c: position of Gaussian
// a: exponential coefficient
// 
// g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
// 
// this version is based on V=V(V1,V2,r1,r2,alpha)
// alpha = (a1 + a2)/(a1 a2)
// 
// dVdr is (1/r)*(dV12/dr)
// dVdV is dV12/dV1 
// dVdalpha is dV12/dalpha
// d2Vdalphadr is (1/r)*d^2V12/dalpha dr
// d2VdVdr is (1/r) d^2V12/dV1 dr
double ogauss_alpha(GaussianVca& g1, GaussianVca& g2, GaussianVca& g12, double& dVdr, double& dVdV, double& sfp);

// an overlap
class GOverlap {
public:
    int level;                      // level (0=root, 1=atoms, 2=2-body, 3=3-body, etc.)
    GaussianVca g;                  // Gaussian representing overlap
    double volume;                  // volume of overlap (also stores Psi1..i in GPU version)
    double dvv1;                    // derivative wrt volume of first atom (also stores F1..i in GPU version)
    double dv1x;                    // derivative wrt position of first atom (also stores P1..i in GPU version) (x)
    double dv1y;                    // derivative wrt position of first atom (also stores P1..i in GPU version) (y)
    double dv1z;                    // derivative wrt position of first atom (also stores P1..i in GPU version) (z)
    double gamma1i;                 // sum gammai for this overlap
    double self_volume;             // self volume accumulator (also stores Psi'1..i in GPU version)
    double sfp;                     // switching function derivatives    
    int atom;                       // the atomic index of the last atom of the overlap list (i, j, k, ..., atom) 
                                    //    = (Parent, atom)
    int parent_index;               // index in tree list of parent overlap
    int children_startindex;        // start index in tree array of children
    int children_count;             // number of children
    void print_overlap(void);
};

// overlap comparison function
bool goverlap_compare(const GOverlap& overlap1, const GOverlap& overlap2);

// A collection of, mainly, recursive routines to constructs and analyze the overlap tree.
// Not meant to be called directly. It is used by GaussVol.
class GOverlap_Tree {
public:
    GOverlap_Tree(int natoms)
    {
        this->natoms = natoms;
    }

    ~GOverlap_Tree(void)
    {
        overlaps.clear();
    }

    //intializes the tree at 1-body level
    int init_overlap_tree(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz, std::vector<double>& radii, std::vector<double>& volumes, std::vector<double>& gammas, std::vector<int>& ishydrogen);

    // adds to the tree the children of overlap identified by "parent_index" in the tree
    int add_children(int parent_index, std::vector<GOverlap>& children_overlaps);

    // scans the siblings of overlap identified by "root_index" to create children overlaps,
    // returns them into the "children_overlaps" buffer: (root) + (atom) -> (root, atom)
    int compute_children(int root_index, std::vector<GOverlap>& children_overlaps);

    //grow the tree with more children starting at the given root slot (recursive)
    int compute_andadd_children_r(int root);

    //compute the tree starting from the 1-body level
    int compute_overlap_tree_r(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz, std::vector<double>& radius, std::vector<double>& volume, std::vector<double>& gamma, std::vector<int>& ishydrogen);

    // compute volumes, energy of the overlap at slot and calls itself recursively to get 
    // the volumes of the children
    int compute_volume_underslot2_r(
        int slot,
        double& psi1i, double& f1i, double& p1ix, double& p1iy, double& p1iz,                            //subtree accumulators for free volume
        double& psip1i, double& fp1i, double& pp1ix, double& pp1iy, double& pp1iz,                       //subtree accumulators for self volume
        double& energy1i, double& fenergy1i, double& penergy1ix, double& penergy1iy, double& penergy1iz, //subtree accumulators for volume-based energy
        std::vector<double>& drx, std::vector<double>& dry, std::vector<double>& drz,                                                                         //gradients of volume-based energy wrt to atomic positions
        std::vector<double>& dv,                                                                         //gradients of volume-based energy wrt to atomic volumes				
        std::vector<double>& free_volume,                                                                //atomic free volumes
        std::vector<double>& self_volume);                                                               //atomic self volumes

    // recursively traverses tree and computes volumes, etc.
    int compute_volume2_r(
        std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz,
        double& volume, double& energy, 
        std::vector<double>& drx, std::vector<double>& dry, std::vector<double>& drz,
        std::vector<double>& dv,
        std::vector<double>& free_volume,
        std::vector<double>& self_volume);

    // rescan the sub-tree to recompute the volumes, does not modify the tree
    int rescan_r(int slot);

    // rescan the tree to recompute the volumes, does not modify the tree
    int rescan_tree_v(
        std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz,
        std::vector<double>& radii,
        std::vector<double>& volumes,
        std::vector<double>& gammas,
        std::vector<int>& ishydrogen);

    // rescan the sub-tree to recompute the gammas, does not modify the volumes nor the tree
    int rescan_gamma_r(int slot);

    // rescan the tree to recompute the gammas only, does not modify volumes and the tree
    int rescan_tree_g(std::vector<double>& gammas);

    //print the contents of the tree
    void print_tree(void);

    //print the contents of the tree (recursive)
    void print_tree_r(int slot);

    //counts number of overlaps under the one given
    int nchildren_under_slot_r(int slot);

    int natoms;
    std::vector<GOverlap> overlaps; //the root is at index 0, atoms are at 1..natoms+1
};

// A class that implements the Gaussian description of an object (molecule) made of a overlapping spheres
class GaussVol {
public: 
    // Creates/Initializes a GaussVol instance
    GaussVol(const int natoms, std::vector<int>& ishydrogen)
    {
        tree = new GOverlap_Tree(natoms);
        this->natoms = natoms;
        this->radii.resize(natoms);
        for(int i=0;i<natoms;i++) radii[i] = 1.;
        this->volumes.resize(natoms);
        for(int i=0;i<natoms;i++) volumes[i] = 0.;
        this->gammas.resize(natoms);
        for(int i=0;i<natoms;i++) gammas[i] = 0.;
        this->ishydrogen = ishydrogen;
    }
    GaussVol(const int natoms, std::vector<double>& radii, std::vector<double>& volumes, std::vector<double>& gammas, std::vector<int>& ishydrogen)
    {
        tree = new GOverlap_Tree(natoms);
        this->natoms = natoms;
        this->radii = radii;
        this->volumes = volumes;
        this->gammas = gammas;
        this->ishydrogen = ishydrogen;
    }
    ~GaussVol(void)
    {
        delete tree;
        radii.clear();
        volumes.clear();
        gammas.clear();
        ishydrogen.clear();
    }

    int setRadii(std::vector<double>& radii)
    {
        if(natoms == radii.size()) {
            this->radii = radii;
            return natoms;
        }
        else {
            throw std::runtime_error("setRadii: number of atoms does not match");
            return -1;
        }
    }

    int setVolumes(std::vector<double>& volumes)
    {
        if(natoms == volumes.size()) {
            this->volumes = volumes;
            return natoms;
        }
        else {
            throw std::runtime_error("setVolumes: number of atoms does not match");
            return -1;
        }
    }

    int setGammas(std::vector<double>& gammas)
    {
        if(natoms == gammas.size()) {
            this->gammas = gammas;
            return natoms;
        }
        else {
            throw std::runtime_error("setGammas: number of atoms does not match");
            return -1;
        }
    }

    //constructs the tree
    void compute_tree(std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz);

    // returns GaussVol volume energy function and forces
    // also returns gradients with respect to atomic volumes and atomic free-volumes and self-volumes
    void compute_volume(
        std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz, 
        double& volume, double& energy, std::vector<double>& forcex, std::vector<double>& forcey, std::vector<double>& forcez,
        std::vector<double>& gradV, std::vector<double>& free_volume,  std::vector<double>& self_volume);

    //rescan the tree after resetting gammas, radii and volumes
    void rescan_tree_volumes(std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz);

    //rescan the tree resetting gammas only with current values
    void rescan_tree_gammas(void);

    // returns number of overlaps for each atom 
    void getstat(std::vector<int>& nov);

    void print_tree(void)
    {
        tree->print_tree();
    }

private:
    GOverlap_Tree *tree;
    int natoms;
    std::vector<double> radii;
    std::vector<double> volumes;
    std::vector<double> gammas;
    std::vector<int> ishydrogen;
};
