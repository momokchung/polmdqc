// Author: Moses KJ Chung
// Year:   2024

#include "hilbert.h"
#include "hilbrt.h"

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  initHilbert  --  initialize the permutation table  //
//                                                     //
/////////////////////////////////////////////////////////

// "initHilbert" initializes the permutation table
//
// The table 'transgc' has 8 x 3 x 8 entries. It contains
// all possible Gray code sequences traveled by the 1st
// order Hilbert curve in 3 dimensions. The first column
// is the Gray code of the entry point of the curve, and
// the second column is the direction (0, 1, or 2, 0 means
// the x-axis) where the exit point of curve lies.
// The table 'tsb1mod3' contains the numbers of trailing
// set '1' bits of the indices from 0 to 7, modulo by '3'.
// The code for generating this table is from:
// https://graphics.stanford.edu/~seander/bithacks.html

void initHilbert(int ndim)
{
    int gc[8],N,mask,travel_bit;
    int c,f,g,k,v;

    for (int i = 0; i < 8; i++) gc[i] = 0;

    N = (ndim == 2) ? 4 : 8;
    mask = (ndim == 2) ? 3 : 7;

    // Generate the Gray code sequence.
    for (int i = 0; i < N; i++) {
        gc[i] = i ^ (i >> 1);
    }

    for (int e = 0; e < N; e++) {
        for (int d = 0; d < ndim; d++) {
            // Calculate the end point (f).
            f = e ^ (1 << d); // Toggle the d-th bit of 'e'.
            // travel_bit = 2**p, the bit we want to travel. 
            travel_bit = e ^ f;
            for (int i = 0; i < N; i++) {
                // // Rotate gc[i] left by (p + 1) % n bits.
                k = gc[i] * (travel_bit * 2);
                g = ((k | (k / N)) & mask);
                // Calculate the permuted Gray code by xor with the start point (e).
                transgc[e][d][i] = (g ^ e);
            }
        } // d
    } // e

    // Count the consecutive '1' bits (trailing) on the right.
    tsb1mod3[0] = 0;
    for (int i = 1; i < N; i++) {
        v = ~i; // Count the 0s.
        v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
        for (c = 0; v; c++) {
            v >>= 1;
        }
        tsb1mod3[i] = c % ndim;
    }
}

inline int splitHilbert(AlfAtom *alfatoms, int size, int gc0, int gc1,
    real xmin, real xmax, real ymin, real ymax, real zmin, real zmax)
{
    int axis,d;
    int i,j;
    real splt;
    AlfAtom swapvert;

    // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which 
    // correspoding to x-, or y- or z-axis.
    axis = (gc0 ^ gc1) >> 1; 

    // Calulate the split position along the axis.
    if (axis == 0) {
        splt = 0.5 *(xmin + xmax);
    }
    else if (axis == 1) {
        splt = 0.5 *(ymin + ymax);
    }
    else { // == 2
        splt = 0.5 *(zmin + zmax);
    }

    // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction
    // of the axis is to the positive of the axis, otherwise, it is -1.
    d = ((gc0 & (1<<axis)) == 0) ? 1 : -1;

    // Partition the vertices into left- and right-arrays such that left points
    // have Hilbert indices lower than the right points.
    i = 0;
    j = size - 1;

    // Partition the vertices into left- and right-arrays.
    if (d > 0) {
        do {
            for (; i < size; i++) {
                if (alfatoms[i].coord[axis] >= splt) break;
            }
            for (; j >= 0; j--) {
                if (alfatoms[j].coord[axis] < splt) break;
            }
            // Is the partition finished?
            if (i == (j + 1)) break;
            // Swap i-th and j-th vertices.
            swapvert = alfatoms[i];
            alfatoms[i] = alfatoms[j];
            alfatoms[j] = swapvert;
            // Continue patitioning the array;
        } while (true);
    }
    else {
        do {
            for (; i < size; i++) {
                if (alfatoms[i].coord[axis] <= splt) break;
            }
            for (; j >= 0; j--) {
                if (alfatoms[j].coord[axis] > splt) break;
            }
            // Is the partition finished?
            if (i == (j + 1)) break;
            // Swap i-th and j-th vertices.
            swapvert = alfatoms[i];
            alfatoms[i] = alfatoms[j];
            alfatoms[j] = swapvert;
            // Continue patitioning the array;
        } while (true);
    }

    return i;
}

/////////////////////////////////////////////////////////////////
//                                                             //
//  sort3DHilbert  --  sort points using the 3d Hilbert curve  //
//                                                             //
/////////////////////////////////////////////////////////////////

// "sort3DHilbert" sorts points using the 3d Hilbert curve

void sort3DHilbert(AlfAtom *alfatoms, int size, int e, int d, real xmin, real ymin, real zmin, real xmax, real ymax, real zmax, int depth)
{
    int N = 3;
    int mask = 7;
    int p[9],e_w,d_w,k,ei,di;
    real x1,x2,y1,y2,z1,z2;

    p[0] = 0;
    p[8] = size;

    // Sort the points according to the 1st order Hilbert curve in 3d.
    p[4] = splitHilbert(&alfatoms[0], p[8], transgc[e][d][3], transgc[e][d][4], xmin, xmax, ymin, ymax, zmin, zmax);
    p[2] = splitHilbert(&alfatoms[0], p[4], transgc[e][d][1], transgc[e][d][2], xmin, xmax, ymin, ymax, zmin, zmax);
    p[1] = splitHilbert(&alfatoms[0], p[2], transgc[e][d][0], transgc[e][d][1], xmin, xmax, ymin, ymax, zmin, zmax);
    p[3] = splitHilbert(&(alfatoms[p[2]]), p[4]-p[2], transgc[e][d][2], transgc[e][d][3], 
        xmin, xmax, ymin, ymax, zmin, zmax) + p[2];
    p[6] = splitHilbert(&(alfatoms[p[4]]), p[8]-p[4], transgc[e][d][5], transgc[e][d][6], 
        xmin, xmax, ymin, ymax, zmin, zmax)+p[4];
    p[5] = splitHilbert(&(alfatoms[p[4]]), p[6]-p[4], transgc[e][d][4], transgc[e][d][5], 
        xmin, xmax, ymin, ymax, zmin, zmax)+p[4];
    p[7] = splitHilbert(&(alfatoms[p[6]]), p[8]-p[6], transgc[e][d][6], transgc[e][d][7], 
        xmin, xmax, ymin, ymax, zmin, zmax)+p[6];

    if (hilbert_order > 0) {
        // A maximum order is prescribed. 
        if ((depth + 1) == hilbert_order) {
            // The maximum prescribed order is reached.
            return;
        }
    }

    // Recursively sort the points in sub-boxes.
    for (int w = 0; w < 8; w++) {
        // w is the local Hilbert index (NOT Gray code).
        // Sort into the sub-box either there are more than 2 points in it, or
        // the prescribed order of the curve is not reached yet.
        //if ((p[w+1] - p[w] > b->hilbert_limit) || (b->hilbert_order > 0)) {
        if ((p[w+1] - p[w]) > hilbert_limit) {
            // Calculcate the start point (ei) of the curve in this sub-box.
            // update e = e ^ (e(w) left_rotate (d+1)).
            if (w == 0) {
                e_w = 0;
            }
            else {
                // calculate e(w) = gc(2 * floor((w - 1) / 2)).
                k = 2 * ((w - 1) / 2); 
                e_w = k ^ (k >> 1); // = gc(k).
            }
            k = e_w;
            e_w = ((k << (d+1)) & mask) | ((k >> (N-d-1)) & mask);
            ei = e ^ e_w;
            // Calulcate the direction (di) of the curve in this sub-box.
            // update d = (d + d(w) + 1) % N
            if (w == 0) {
                d_w = 0;
            }
            else {
                d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
            }
            di = (d + d_w + 1) % N;
            // Calculate the bounding box of the sub-box.
            if (transgc[e][d][w] & 1) { // x-axis
                x1 = 0.5 * (xmin + xmax);
                x2 = xmax;
            }
            else {
                x1 = xmin;
                x2 = 0.5 * (xmin + xmax);
            }
            if (transgc[e][d][w] & 2) { // y-axis
                y1 = 0.5 * (ymin + ymax);
                y2 = ymax;
            }
            else {
                y1 = ymin;
                y2 = 0.5 * (ymin + ymax);
            }
            if (transgc[e][d][w] & 4) { // z-axis
                z1 = 0.5 * (zmin + zmax);
                z2 = zmax;
            }
            else {
                z1 = zmin;
                z2 = 0.5 * (zmin + zmax);
            }
            sort3DHilbert(&(alfatoms[p[w]]), p[w+1] - p[w], ei, di, x1, y1, z1, x2, y2, z2, depth+1);
        } // if (p[w+1] - p[w] > 1)
    } // w
}

//////////////////////////////////////////////////////////////////
//                                                              //
//  brioHilbert  --  sort points using brio / 3d Hilbert curve  //
//                                                              //
//////////////////////////////////////////////////////////////////

// "brioHilbert" sorts points using brio / 3d Hilbert curve

void brioHilbert(AlfAtom *alfatoms, int size, real xmin, real ymin, real zmin, real xmax, real ymax, real zmax, int depth)
{
    int middle = 0;

    if(size >= brio_threshold) {
        depth++;
        middle = size*brio_ratio;
        brioHilbert(alfatoms, middle, xmin, ymin, zmin, xmax, ymax, zmax, depth);
    }

    sort3DHilbert(&(alfatoms[middle]), size - middle, 0, 0, xmin, ymin, zmin, xmax, ymax, zmax, 0);
}
}
