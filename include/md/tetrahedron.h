/* ===========================================================================================
*
*	For each tetrahedron, store:
*
*	- the index of its four vertices
*	- pointers to its neighbours (4 maximum).
*		neighbor(i) is the tetrahedron that shares
*		all vertices of the tetrahedron considered, except i
*		(0 if the corresponding face is on the convex hull)
*
*
*       Tetrahedron:    define two arrays of integer*1:
*                       tetra_info stores:
*                               bit 0: orientation
*                               bit 1: status
*                               bit 2-5: surface info 
*					(one for each face of the tetrahedron)
*					surface info is a tag on the
*					face of the tetrahedron considered.
*					This tag can be set to represent
*					the convex hull of the molecule
*					for example, or to indicate if the
*					face belongs to a restricted Delaunay
*					(such as the one used in a skin
*					surface)
*                       tetra_nindex:
*				if tetra is (a,b,c,d), each of its
*				face is shared with another tetrahedron.
*				For example, face (b,c,d) is also a
*				face of tetrahedron (b,c,d,e). 
*				tetra_nindex gives the index of point e
*				in the list (b,c,d,e)
*				This index can take value 0, 1, 2, 3
*
*
 ==============================================================================================*/

#pragma once

#include "vertex.h"

class Tetrahedron {
public:
	int Vertices[4]; 		// the four vertices defining the tetrahedron
	int Neighbours[4];		// the four neighbours (-1 if no neighbour)
	std::bitset<8> info;		// (see above)
	int info_edge[6];
	short nindex[4];

	Tetrahedron();
	~Tetrahedron();

	void init();

};
