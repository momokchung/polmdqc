/* ====================================================================
   delcx
 
	This class computes the regular triangulation of a set
	of N weighted points in 3D, using the incremental flipping
	algorithm of Edelsbrunner.

	This implementation is based on the algorithm published in
	H. Edelsbrunner and N.R. Shah, Algorithmica (1996) 15: 223-241

	1) Algorithm:
	*************

	Briefly, the algorithm works as follows:

	- first, a large tetrahedron initialises the program. All
	four vertices of this tetrahedron are set at "infinite"

	- All N points are added one by one.

	- For each point:

		- localize the tetrahedron in the current regular
		triangulation that contains this point

		- test if the point is redundant; if yes, remove

		- If the point is not redundant, insert in the
		tetrahedron : this is a "1-4" flip

		- collect all "link facets" (i.e. all triangles
		in tetrahedron containing the new point, that face
		this new point) that are not regular.

		- for each non-regular link facet, check if it
		is "flippable". If yes, perform a "2-3", "3-2"
		or "1-4" flip. Add new link facets in the list,
		if needed.

		- when link facet list if empty, move to next 
		point

	- Remove "infinite" tetrahedron, i.e. tetrahedron with
	one vertice at "infinite"

	- collect all remaining tetrahedron, define convex hull,
	and exit.

	2) Data structure:
	******************

	I maintain a minimal data structure that includes only 
	the tetrahedrons of the triangulation (triangles,
	edges and vertices are implicit).

	For each tetrahedron, I store:

	- the index of its four vertices
	- pointers to its neighbours (4 maximum).
		neighbor(i) is the tetrahedron that shares
		all vertices of the tetrahedron considered, except i
		(0 if the corresponding face is on the convex hull)
	- its status: 1 "active" (i.e. part of the triangulation), 0 inactive
	- its orientation


	2) number representation:
	**************************

	I use double precision floating points. However, if one of
	the geometricaly test becomes "imprecise", I switch to
	arbitrary precision arithmetics (using the gmp package).

	All predicates have therefore a floating point filter

 ==================================================================== */

#pragma once
#include "vertex.h"
#include "vector.h"
#include "tetrahedron.h"
#include "sosgmp.h"
#include "sorttools.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>
#include <stack>

/*********************************************************************************
  Delcx class
 *********************************************************************************/

class DELCX {

public:
	// Setup calculations
	void setup(int npoints, double *coord, double *radii, double *coefS, double *coefV,
	double *coefM, double *coefG, std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// Compute 3D weighted Delaunay triangulation
	void regular3D(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// generate list of edges
	void delaunayEdges(std::vector<Tetrahedron>& tetra, std::vector<std::pair<int, int> >& edges);

private:

	// locate tetrahedron in which point is inserted
	void locate_jw(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, 
			int ipoint, int *tetra_loc, int *iredundant);

	// go over full link_facet
	void flip(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// Check if a point is inside a tetrahedron
	void inside_tetra(std::vector<Vertex>& vertices, int p, int a, int b, int c, int d, 
		int iorient, bool *is_in, bool *redundant, int *ifail);
 
	// Check if a facet connects two tetrehedra that are convex
	void regular_convex(std::vector<Vertex>& vertices, int a, int b, int c, 
		int p, int o, int itest_abcp, bool *regular, bool *convex, bool *test_abpo, 
		bool *test_bcpo, bool *test_capo);

	// sign associated with the missing inf point
	void missinf_sign(int i, int j, int k, int *l, int *sign);

	// flip 1-4: from 1 to 4 tetrahedra
	void flip_1_4(std::vector<Tetrahedron>& tetra, int ipoint, int itetra, 
		int *tetra_last);

	// flip 4-1: from 4 to 1 tetrahedron
	void flip_4_1(std::vector<Vertex>& vertice, std::vector<Tetrahedron>& tetra, 
	int itetra, int jtetra, int ktetra, int ltetra, int *vertices,
	int idp, int jdp, int kdp, int ldp, bool test_acpo, int *ierr, int *tetra_last);

	// flip 2-3: from 2 to 3 tetrahedra
	void flip_2_3(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
	int *vertices, int *facei, int *facej, bool test_abpo, bool test_bcpo, 
	bool test_capo, int *ierr, int *tetra_last);

	// flip 3-2: from 3 to 2 tetrahedra
	void flip_3_2(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int ktetra, int *vertices, int *edgei, int *edgej, int *edgek,
		bool test_bcpo, bool test_acpo, int *ierr, int *tetra_last);

	// Define facet between two tetrahedron
	void define_facet(std::vector<Tetrahedron>& tetra, int itetra, int jtetra, 
		int idx_o, int *facei, int *facej);

	// info about tetrahedron
	void find_tetra(std::vector<Tetrahedron>& tetra, int itetra, int idx_c, 
	int a, int b, int o, int *ifind, int *tetra_loc, int *idx_a, int *idx_b);

	// reorder vertices of tetrahedra in increasing order
	void reorder_tetra(std::vector<Tetrahedron>& tetra);

	// remove "infinite" tetrahedron
	void remove_inf(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

	// mark tetrahedron
	void mark_zero(std::vector<Tetrahedron>& tetra, int itetra, int ivertex);

	// Computes the volume of a tetrahedron
	double tetra_vol(double *a, double *b, double *c, double *d);

	// remove flat tetrahedra
	void peel(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra, int *flag);

	// find edge in tetra defined by 2 vertices
	int findEdge(Tetrahedron t, int i1, int i2);

	// Add bogus points as needed so that we have at least 4 points
	void addBogus(int npoints, double *coord, double *radii, double *bcoord, double *brad); 

protected:

	std::queue<std::pair<int, int> > link_facet;
	std::queue<std::pair<int, int> > link_index;
	std::stack<int> free;
	std::vector<int> kill;

	double eps = 1.e-3;
	double eps_vol = 1e-3;

	int inf4_1[4] = {1, 1, 0, 0};
	int sign4_1[4] = {-1, 1, 1, -1};
	int inf4_2[4][4] = {
		{ -1, 1, 2, 2},
		{ 1, -1, 2, 2},
		{ 2, 2, -1, 0},
		{ 2, 2, 0, -1}
	};
	int sign4_2[4][4] = {
		{ 0, 1, -1, 1},
		{ -1, 0, 1, -1},
		{ 1, -1, 0, 1},
		{ -1, 1, -1, 0}
	};
	int sign4_3[4] = {-1, 1, -1, 1};
	int inf5_2[4][4] = {
		{ -1, 1, 0, 0},
		{ 1, -1, 0, 0},
		{ 0, 0, -1, 0},
		{ 0, 0, 0, -1}
	};
	int sign5_2[4][4] = {
		{ 0, -1, -1, 1},
		{ 1, 0, -1, 1},
		{ 1, 1, 0, 1},
		{ -1, -1, -1, 0}
	};
	int inf5_3[4] = {0, 0, 2, 2};
	int sign5_3[4] = {1, 1, -1, 1};
	int order1[4][3] = {
		{ 2, 1, 3},
		{ 0, 2, 3},
		{ 1, 0, 3},
		{ 0, 1, 2}
	};

	int ord_rc[3][3] = {
		{0, 1, 2},
		{2, 0, 1},
		{1, 2, 0},
	};

	int order2[6][2] = {
		{ 2, 3},
		{ 3, 1},
		{ 1, 2},
		{ 0, 3},
		{ 2, 0},
		{ 0, 1}
	};
	int order3[6][2] = {
		{ 0, 1},
		{ 0, 2},
		{ 0, 3},
		{ 1, 2},
		{ 1, 3},
		{ 2, 3}
	};
	int idxList[4][3] = {
		{ 0, 0, 0},
		{ 0, 1, 1},
		{ 1, 1, 2},
		{ 2, 2, 2}
	};
	int table32[3][3] = {
		{ 0, 1, 2},
		{ 0, 2, 1},
		{ 2, 0, 1}
	};
	int table32_2[3][2] = {
		{ 0, 1},
		{ 0, 2},
		{ 1, 2}
	};
	int table41[3][3] = {
		{ 1, 0, 2},
		{ 0, 1, 2},
		{ 0, 2, 1}
	};
	int table41_2[3][2] = {
		{ 0, 0},
		{ 1, 0},
		{ 1, 1}
	};
	int order[3][2] = {
		{ 1, 2},
		{ 2, 0},
		{ 0, 1}
	};
	int other[4][3] = {
		{ 1, 2, 3},
		{ 0, 2, 3},
		{ 0, 1, 3},
		{ 0, 1, 2}
	};
	int other2[4][4][2] = {
		{
			{ -1, -1},
			{ 2, 3},
			{ 1, 3},
			{ 1, 2}
		},
		{
			{ 2, 3},
			{ -1, -1},
			{ 0, 3},
			{ 0, 2}
		},
		{
			{ 1, 3},
			{ 0, 3},
			{ -1, -1},
			{ 0, 1}
		},
		{
			{ 1, 2},
			{ 0, 2},
			{ 0, 1},
			{ -1, -1}
		},
	};
};
