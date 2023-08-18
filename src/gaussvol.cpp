#include "gaussvol.h"
#include "mathConst.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

//counts overlaps
static int _nov_ = 0;

// overlap volume switching function + 1st derivative
double pol_switchfunc(double gvol, double volmina, double volminb, double& sp)
{
    double swf = 0.;
    double swfp = 1.;
    double swd, swu, swu2, swu3, s;
    if(gvol > volminb) {
        swf = 1.;
        swfp = 0.;
    }
    else if (gvol < volmina) {
        swf = 0.;
        swfp = 0.;
    }
    swd = 1./(volminb - volmina);
    swu = (gvol - volmina)*swd;
    swu2 = swu*swu;
    swu3 = swu*swu2;
    s = swf + swfp*swu3*(10.-15.*swu+6.*swu2);
    sp = swfp*swd*30.*swu2*(1. - 2.*swu + swu2);

    //turn off switching function
    //*sp = 0.0;
    //s = 1.0;
    return s;
}

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
double ogauss_alpha(GaussianVca& g1, GaussianVca& g2, GaussianVca& g12, double& dVdr, double& dVdV, double& sfp)
{
    double d2, dx, dy, dz;
    double c1x = g1.cx;
    double c1y = g1.cy;
    double c1z = g1.cz;
    double c2x = g2.cx;
    double c2y = g2.cy;
    double c2z = g2.cz;
    double distx;
    double disty;
    double distz;
    double deltai, gvol, p12, a12;
    double s, sp, df, dgvol, dgvolv, ef, dgvola2, dgvola1, dgalpha, dgalpha2, dgvolvdr;

    distx = c2x - c1x;
    disty = c2y - c1y;
    distz = c2z - c1z;
    d2 = distx*distx + disty*disty + distz*distz;

    a12 = g1.a + g2.a;
    deltai = 1./a12;
    df = (g1.a)*(g2.a)*deltai; // 1/alpha

    ef = std::exp(-df*d2);
    gvol = ((g1.v * g2.v)/std::pow(pi/df,1.5))*ef;
    dgvol = -2.*df*gvol; // (1/r)*(dV/dr) w/o switching function
    dgvolv = g1.v > 0. ? gvol/g1.v : 0.;     // (dV/dV1)  w/o switching function

    // parameters for overlap gaussian. Note that c1 and c2 are Vec3's and the "*" operator wants 
    // the vector first and scalar second vector2 = vector1 * scalar
    g12.cx = ((c1x * g1.a) + (c2x * g2.a)) * deltai;
    g12.cy = ((c1y * g1.a) + (c2y * g2.a)) * deltai;
    g12.cz = ((c1z * g1.a) + (c2z * g2.a)) * deltai;
    g12.a = a12;
    g12.v = gvol;

    // switching function
    s = pol_switchfunc(gvol, volmina, volminb, sp);
    sfp = sp*gvol+s;
    dVdr = dgvol;
    dVdV = dgvolv;

    return s*gvol;
}

// overlap comparison function
bool goverlap_compare(const GOverlap& overlap1, const GOverlap& overlap2)
{
    // order by volume, larger first
    return overlap1.volume > overlap2.volume;
}

int GOverlap_Tree::init_overlap_tree(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz, std::vector<double>& radius, std::vector<double>& volume, std::vector<double>& gamma, std::vector<int>& ishydrogen)
{
    GOverlap overlap;

    // reset tree
    overlaps.clear();

    // slot 0 contains the master tree information, children = all of the atoms
    overlap.level = 0;
    overlap.volume = 0;
    overlap.dv1x = 0.;
    overlap.dv1y = 0.;
    overlap.dv1z = 0.;
    overlap.dvv1 = 0.;
    overlap.self_volume = 0;
    overlap.sfp = 1.;
    overlap.gamma1i = 0.;
    overlap.parent_index = -1;
    overlap.atom = -1;
    overlap.children_startindex = 1;
    overlap.children_count = natoms;

    overlaps.push_back(overlap);

    // list of atoms start at slot #1
    for(int iat = 0; iat < natoms; iat++){
        double a = kfc/(radius[iat]*radius[iat]);
        double vol = ishydrogen[iat] > 0 ? 0. : volume[iat];
        overlap.level = 1;
        overlap.g.v = vol;
        overlap.g.a = a;
        overlap.g.cx = posx[iat];
        overlap.g.cy = posy[iat];
        overlap.g.cz = posz[iat];
        overlap.volume = vol;
        overlap.dv1x = 0.;
        overlap.dv1y = 0.;
        overlap.dv1z = 0.;
        overlap.dvv1 = 1.; //dVi/dVi
        overlap.self_volume = 0.;
        overlap.sfp = 1.;
        overlap.gamma1i = gamma[iat];// gamma[iat]/SA_DR;
        overlap.parent_index = 0;
        overlap.atom = iat; 
        overlap.children_startindex = -1;
        overlap.children_count = -1;
        overlaps.push_back(overlap);
    }

    return 1;
}

// adds to the tree the children of overlap identified by "parent_index" in the tree
int GOverlap_Tree::add_children(int parent_index, std::vector<GOverlap>& children_overlaps)
{
    int i, ip, slot;

    // adds children starting at the last slot
    int start_index = overlaps.size();

    int noverlaps = children_overlaps.size();

    // retrieves address of root overlap
    GOverlap* root = &(overlaps[parent_index]);

    // registers list of children
    root->children_startindex = start_index;
    root->children_count = noverlaps;

    // sort neighbors by overlap volume
    //if(root->level == 1){
    std::sort(children_overlaps.begin(), children_overlaps.end(), goverlap_compare);
        //}

    int root_level = root->level;

    // now copies the children overlaps from temp buffer
    for(int ip = 0; ip < noverlaps; ip++){
        children_overlaps[ip].level = root_level + 1;
        // connect overlap to parent
        children_overlaps[ip].parent_index = parent_index;
        // reset its children indexes 
        children_overlaps[ip].children_startindex = -1;
        children_overlaps[ip].children_count = -1;
        // add to tree
        // note that the 'root' pointer may be invalidated by the push back below
        overlaps.push_back(children_overlaps[ip]);
    }

    _nov_ += noverlaps;

    return start_index;
}


// scans the siblings of overlap identified by "root_index" to create children overlaps,
// returns them into the "children_overlaps" buffer: (root) + (atom) -> (root, atom)
int GOverlap_Tree::compute_children(int root_index, std::vector<GOverlap>& children_overlaps)
{
    int parent_index;
    int sibling_start, sibling_count;
    int j;

    // reset output buffer
    children_overlaps.clear();

    // retrieves overlap
    GOverlap& root = overlaps[root_index];
  
    // retrieves parent overlap
    parent_index = root.parent_index;
    if(parent_index < 0) return 1; // master root? can't do compute_children() on master root
    if(root.level >= max_order) return 1;
    GOverlap& parent = overlaps[parent_index];

    // retrieves start index and count of siblings
    sibling_start = parent.children_startindex;
    sibling_count = parent.children_count;
    if(sibling_start < 0 || sibling_count < 0) return -1; // parent is not initialized?
    if(root_index < sibling_start && root_index > sibling_start + sibling_count -1 ) return -1; // this overlap somehow is not the child of registered parent

    // now loops over "younger" siblings (i<j loop) to compute new overlaps
    for(int slotj = root_index+1; slotj < sibling_start+sibling_count; slotj++){
        GaussianVca g12;
        GOverlap& sibling = overlaps[slotj];
        double gvol, dVdr,dVdV, sfp;

        // atomic gaussian of last atom of sibling
        int atom2 = sibling.atom;
        GaussianVca& g1 = root.g;
        GaussianVca& g2 = overlaps[atom2+1].g; //atoms are stored in the tree at indexes 1...N
        gvol = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp);

        // create child if overlap volume is not zero
        if(gvol > min_gvol){
            GOverlap ov;
            ov.g = g12;
            ov.volume = gvol;
            ov.self_volume = 0;
            ov.atom = atom2;
            // dv1 is the gradient of V(123..)n with respect to the position of 1
            ov.dv1x = (g2.cx - g1.cx) * (-dVdr);
            ov.dv1y = (g2.cy - g1.cy) * (-dVdr);
            ov.dv1z = (g2.cz - g1.cz) * (-dVdr);
            //dvv1 is the derivative of V(123...)n with respect to V(123...)
            ov.dvv1 = dVdV;
            ov.sfp = sfp;
            ov.gamma1i = root.gamma1i + overlaps[atom2+1].gamma1i;
            children_overlaps.push_back(ov);
        }
    }

    return 1;
}


// rescan the sub-tree to recompute the volumes, does not modify the tree
int GOverlap_Tree::rescan_r(int slot)
{
    int parent_index;
    int sibling_start, sibling_count;

    // this overlap
    GOverlap& ov = overlaps[slot];

    // recompute its own overlap by merging parent and last atom
    parent_index = ov.parent_index;
    if (parent_index > 0) {
        GaussianVca g12;
        double dVdr,dVdV, dVdalpha, d2Vdalphadr, d2VdVdr, sfp;
        int atom = ov.atom;
        GOverlap& parent  = overlaps[parent_index];
        GaussianVca& g1 = parent.g;
        GaussianVca& g2 = overlaps[atom+1].g; //atoms are stored in the tree at indexes 1...N
        double gvol = ogauss_alpha(g1,g2, g12,dVdr,dVdV,sfp);
        ov.g = g12;
        ov.volume = gvol;
        // dv1 is the gradient of V(123..)n with respect to the position of 1
        ov.dv1x = (g2.cx - g1.cx) * (-dVdr);
        ov.dv1y = (g2.cy - g1.cy) * (-dVdr);
        ov.dv1z = (g2.cz - g1.cz) * (-dVdr);
        //dvv1 is the derivative of V(123...)n with respect to V(123...)
        ov.dvv1 = dVdV;
        ov.sfp = sfp;
        ov.gamma1i = parent.gamma1i + overlaps[atom+1].gamma1i;
    }

    // calls itself recursively on the children
    for(int slot_child = ov.children_startindex; slot_child < ov.children_startindex+ov.children_count; slot_child++) {
        rescan_r(slot_child);
    }

    return 1;
}

// rescan the tree to recompute the volumes, does not modify the tree
int GOverlap_Tree::rescan_tree_v(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz, std::vector<double>& radius, std::vector<double>& volume, std::vector<double>& gamma, std::vector<int>& ishydrogen)
{
    int slot;

    slot = 0;
    GOverlap* ov = &(overlaps[slot]);
    ov->level = 0;
    ov->volume = 0;
    ov->dv1x = 0.;
    ov->dv1y = 0.;
    ov->dv1z = 0.;
    ov->dvv1 = 0.;
    ov->self_volume = 0;
    ov->sfp = 1.;
    ov->gamma1i = 0.;

    slot = 1;
    for(int iat = 0; iat < natoms; iat++, slot++){
        double a = kfc/(radius[iat]*radius[iat]);
        double vol = ishydrogen[iat] > 0 ? 0. : volume[iat];
        ov = &(overlaps[slot]);
        ov->level = 1;
        ov->g.v = vol;
        ov->g.a = a;
        ov->g.cx = posx[iat];
        ov->g.cy = posy[iat];
        ov->g.cz = posz[iat];
        ov->volume = vol;
        ov->dv1x = 0.;
        ov->dv1y = 0.;
        ov->dv1z = 0.;
        ov->dvv1 = 1.; //dVi/dVi
        ov->self_volume = 0.;
        ov->sfp = 1.;
        ov->gamma1i = gamma[iat]; // gamma[iat]/SA_DR;
    }

    rescan_r(0);
    return 1;
}

// rescan the sub-tree to recompute the gammas, does not modify the volumes nor the tree
int GOverlap_Tree::rescan_gamma_r(int slot)
{
    int parent_index;
    int sibling_start, sibling_count;

    // this overlap
    GOverlap& ov = overlaps[slot];

    // recompute its own overlap by merging parent and last atom
    parent_index = ov.parent_index;
    if(parent_index > 0) {
        int atom = ov.atom;
        GOverlap& parent  = overlaps[parent_index];
        ov.gamma1i = parent.gamma1i + overlaps[atom+1].gamma1i;
    }
    
    /* calls itself recursively on the children */
    for(int slot_child = ov.children_startindex; slot_child < ov.children_startindex+ov.children_count; slot_child++) {
        rescan_gamma_r(slot_child);
    }

    return 1;
}

// rescan the tree to recompute the gammas only, does not modify volumes and the tree
int GOverlap_Tree::rescan_tree_g(std::vector<double>& gamma)
{
    int slot;

    slot = 0;
    GOverlap* ov = &(overlaps[slot]);
    ov->gamma1i = 0.;

    slot = 1;
    for(int iat = 0; iat < natoms; iat++, slot++){
        ov = &(overlaps[slot]);
        ov->gamma1i = gamma[iat];
    }

    rescan_gamma_r(0);
    return 1;
}



int GOverlap_Tree::compute_andadd_children_r(int root)
{
    std::vector<GOverlap> children_overlaps;
    compute_children(root, children_overlaps);
    int noverlaps = children_overlaps.size();
    if(noverlaps > 0) {
        int start_slot = add_children(root, children_overlaps);
        for (int ichild=start_slot; ichild < start_slot+noverlaps ; ichild++) {
            compute_andadd_children_r(ichild);
        }
    }
    return 1;
}

int GOverlap_Tree::compute_overlap_tree_r(std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz, std::vector<double>& radius, std::vector<double>& volume, std::vector<double>& gamma, std::vector<int>& ishydrogen)
{
    init_overlap_tree(posx, posy, posz, radius, volume, gamma, ishydrogen);
    for(int slot = 1; slot <= natoms ; slot++){
        compute_andadd_children_r(slot);
    }
    return 1;
}

// compute volumes, energy of this volume and calls itself to get the volumes of the children
int GOverlap_Tree::compute_volume_underslot2_r(
    int slot,
    double& psi1i, double& f1i, double& p1ix, double& p1iy, double& p1iz,                            //subtree accumulators for free volume
    double& psip1i, double& fp1i, double& pp1ix, double& pp1iy, double& pp1iz,                       //subtree accumulators for self volume
    double& energy1i, double& fenergy1i, double& penergy1ix, double& penergy1iy, double& penergy1iz, //subtree accumulators for volume-based energy
    std::vector<double>& drx, std::vector<double>& dry, std::vector<double>& drz,                    //gradients of volume-based energy wrt to atom positions
    std::vector<double>& dv,                                                                         //gradients of volume-based energy wrt to atomic volumes
    std::vector<double>& free_volume,                                                                //atomic free volumes
    std::vector<double>& self_volume)                                                                //atomic self volumes
{
    GOverlap& ov = overlaps[slot];
    double cf = ov.level % 2 == 0 ? -1.0 : 1.0;
    double volcoeff  = ov.level > 0 ? cf : 0;
    double volcoeffp = ov.level > 0 ? volcoeff/(double)ov.level : 0;

    int atom = ov.atom;
    double ai = overlaps[atom+1].g.a;
    double a1i = ov.g.a;
    double a1 = a1i - ai;

    int i,j;
    double c1, c1p, c2;

    psi1i = volcoeff*ov.volume; //for free volumes
    f1i = volcoeff*ov.sfp ;
    p1ix = 0.;
    p1iy = 0.;
    p1iz = 0.;

    psip1i = volcoeffp*ov.volume; //for self volumes
    fp1i = volcoeffp*ov.sfp;
    pp1ix = 0.;
    pp1iy = 0.;
    pp1iz = 0.;

    energy1i = volcoeffp*ov.gamma1i*ov.volume; //EV energy
    fenergy1i = volcoeffp*ov.sfp*ov.gamma1i;
    penergy1ix = 0.;
    penergy1iy = 0.;
    penergy1iz = 0.;

    if(ov.children_startindex >= 0) {
        int sloti;
        for(int sloti=ov.children_startindex; sloti < ov.children_startindex+ov.children_count; sloti++) {
            double psi1it, f1it, p1itx, p1ity, p1itz;
            double psip1it, fp1it, pp1itx, pp1ity, pp1itz;
            double energy1it, fenergy1it, penergy1itx, penergy1ity, penergy1itz;
            compute_volume_underslot2_r(sloti,
                        psi1it, f1it, p1itx, p1ity, p1itz,
                        psip1it, fp1it, pp1itx, pp1ity, pp1itz,
                        energy1it, fenergy1it, penergy1itx, penergy1ity, penergy1itz,
                        drx, dry, drz, dv, free_volume, self_volume);

            psi1i += psi1it;
            f1i += f1it;
            p1ix += p1itx;
            p1iy += p1ity;
            p1iz += p1itz;

            psip1i += psip1it;
            fp1i += fp1it;
            pp1ix += pp1itx;
            pp1iy += pp1ity;
            pp1iz += pp1itz;

            energy1i += energy1it;
            fenergy1i += fenergy1it;
            penergy1ix += penergy1itx;
            penergy1iy += penergy1ity;
            penergy1iz += penergy1itz;
        }
    }

    if(ov.level > 0){
        //contributions to free and self volume of last atom
        free_volume[atom] += psi1i;
        self_volume[atom] += psip1i;

        //contributions to energy gradients
        c2 = ai/a1i;
        drx[atom] += (-ov.dv1x) * fenergy1i + penergy1ix * c2;
        dry[atom] += (-ov.dv1y) * fenergy1i + penergy1iy * c2;
        drz[atom] += (-ov.dv1z) * fenergy1i + penergy1iz * c2;
        //ov.g.v is the unswitched volume
        dv[atom] += ov.g.v * fenergy1i; //will be divided by Vatom later 
        
        //update subtree P1..i's for parent
        c2 = a1/a1i;
        p1ix = (ov.dv1x) * f1i + p1ix * c2;
        p1iy = (ov.dv1x) * f1i + p1iy * c2;
        p1iz = (ov.dv1x) * f1i + p1iz * c2;
        pp1ix = (ov.dv1x) * fp1i + pp1ix * c2;
        pp1iy = (ov.dv1y) * fp1i + pp1iy * c2;
        pp1iz = (ov.dv1z) * fp1i + pp1iz * c2;
        penergy1ix = (ov.dv1x) * fenergy1i + penergy1ix * c2;
        penergy1iy = (ov.dv1y) * fenergy1i + penergy1iy * c2;
        penergy1iz = (ov.dv1z) * fenergy1i + penergy1iz * c2;
        //update subtree F1..i's for parent
        f1i = ov.dvv1 * f1i;
        fp1i = ov.dvv1 * fp1i;
        fenergy1i = ov.dvv1 * fenergy1i;
    }
    return 1;
}

/* traverses tree and computes volumes, etc. */
int GOverlap_Tree::compute_volume2_r(
    std::vector<double>& posx, std::vector<double>& posy, std::vector<double>& posz,
	double& volume, double& energy, 
	std::vector<double>& drx, std::vector<double>& dry, std::vector<double>& drz,
	std::vector<double>& dv,
	std::vector<double>& free_volume,
	std::vector<double>& self_volume)
{ 
    int slot = 0;
    int i,j;
    double psi1i, f1i, p1ix, p1iy, p1iz;                            //subtree accumulators for (free) volume
    double psip1i, fp1i, pp1ix, pp1iy, pp1iz;                       //subtree accumulators for self volume
    double energy1i, fenergy1i, penergy1ix, penergy1iy, penergy1iz; //subtree accumulators for volume-based energy

    // reset volumes, gradients
    for(int i = 0 ; i < drx.size(); ++i) drx[i] = 0.;
    for(int i = 0 ; i < dry.size(); ++i) dry[i] = 0.;
    for(int i = 0 ; i < drz.size(); ++i) drz[i] = 0.;
    for(int i = 0 ; i < dv.size(); ++i) dv[i] = 0.;
    for(int i = 0 ; i < free_volume.size(); ++i) free_volume[i] = 0;
    for(int i = 0 ; i < self_volume.size(); ++i) self_volume[i] = 0;

    compute_volume_underslot2_r(
        0,
        psi1i, f1i, p1ix, p1iy, p1iz,
        psip1i, fp1i, pp1ix, pp1iy, pp1iz,
        energy1i, fenergy1i, penergy1ix, penergy1iy, penergy1iz,
        drx, dry, drz, dv, free_volume, self_volume);

    volume = psi1i;
    energy = energy1i;
    return 1;
}

#ifdef NOTNOW
// print overlaps up to 2-body
static void print_flat_tree_2body(GOverlap_Tree& tree)
{
    int end = 0;
    // the end of the 2-body must be given by the last child of the last atom with children
    for(int slot = 1; slot <= tree.natoms; slot++) {
        if(tree.overlaps[slot].children_startindex > 0) end = tree.overlaps[slot].children_startindex;
    }
    //now print
    for(int slot = 0; slot <= end ; slot++) {
        GOverlap* ov = &(tree.overlaps[slot]);
        std::cout << slot << " " << ov->volume << " " << ov->children_startindex << " " << ov->children_count << std::endl;
    }
}
#endif

#ifdef NOTNOW
static void test_gaussian(GOverlap_Tree& tree)
{
    // test ogauss derivatives
    GaussianVca g1, g2, g12;
    int iter, niter = 1000;
    double d, dx, gvol, dVdr, dVdV, gvol_old, sfp;
    double gradx, grady, gradz;
    double distx, disty, distz;
    double dist2;

    g1 = tree.overlaps[1].g;
    g2 = tree.overlaps[2].g;
    gvol = gvol_old = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp);
    distx = g2.cx - g1.cx;
    disty = g2.cy - g1.cy;
    distz = g2.cz - g1.cz;
    dist2 = distx*distx + disty*disty + distz*distz;
    gradx = distx * (-dVdr) * sfp; //gradient with respect to position of g2
    grady = disty * (-dVdr) * sfp;
    gradz = distz * (-dVdr) * sfp;
    d = std::sqrt(dist2);
    dx = 0.01;
    for(iter = 0; iter < niter; iter++){
        g1.cz += dx;
        gvol = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp);
        distx = g2.cx - g1.cx;
        disty = g2.cy - g1.cy;
        distz = g2.cz - g1.cz;
        dist2 = distx*distx + disty*disty + distz*distz;
        gradx = distx * (-dVdr) * sfp; //gradient with respect to position of g2
        grady = disty * (-dVdr) * sfp;
        gradz = distz * (-dVdr) * sfp;
        d = std::sqrt(dist2);
        std::cout << d << " " << gvol << " " << gvol-gvol_old << " " << gradz*dx << std::endl;
        gvol_old = gvol;
    }
}
#endif

void GOverlap::print_overlap(void)
{
    std::cout << std::setprecision(4) << std::setw(7) << level << " " << std::setw(7)  << atom << " " << std::setw(7)  << parent_index << " " <<  std::setw(7) << children_startindex << " " << std::setw(7) << children_count << " " << std::setw(10) << self_volume << " " << std::setw(10) << volume << " " << std::setw(10) << gamma1i << " " << std::setw(10) << g.a << " " << std::setw(10) << g.cx << " " <<  std::setw(10) << g.cy << " " <<  std::setw(10) << g.cz << " " <<  std::setw(10) << dv1x << " " << std::setw(10) << dv1y << " " << std::setw(10) << dv1z << " " << std::setw(10) << sfp << std::endl;
}

void GOverlap_Tree::print_tree_r(int slot)
{
    GOverlap& ov = overlaps[slot];
    std::cout << "tg: " << std::setw(6) << slot << " ";
    ov.print_overlap();
    for(int i = ov.children_startindex; i < ov.children_startindex+ov.children_count; i++){
        print_tree_r(i);
    }
}

void GOverlap_Tree::print_tree(void)
{
    std::cout << "slot level LastAtom parent ChStart ChCount SelfV V gamma a x y z dedx dedy dedz sfp" << std::endl;
    for(int i = 1; i<= natoms; i++){
        print_tree_r(i);
    }
}

void GaussVol::compute_tree(std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz)
{
    tree->compute_overlap_tree_r(positionsx, positionsy, positionsz, radii, volumes, gammas, ishydrogen);
}

void GaussVol::compute_volume(
    std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz,
    double& volume,
    double& energy,
    std::vector<double>& forcex, std::vector<double>& forcey, std::vector<double>& forcez,
    std::vector<double>& gradV,
    std::vector<double>& free_volume,  std::vector<double>& self_volume)
{
    tree->compute_volume2_r(
        positionsx, positionsy, positionsz,
        volume, energy, 
        forcex, forcey, forcez,
        gradV,
        free_volume, self_volume); 
    for(int i = 0; i < natoms; ++i) {
        forcex[i] = -forcex[i]; //transform gradient to force
        forcey[i] = -forcey[i];
        forcez[i] = -forcez[i];
    }
    for(int i = 0; i < natoms; ++i) {
        if(volumes[i] > 0) gradV[i] = gradV[i]/volumes[i];
    }
}

//rescan to compute a subset of overlap volumes with radii smaller than ones used to
//set up the tree with compute_tree()
void GaussVol::rescan_tree_volumes(std::vector<double>& positionsx, std::vector<double>& positionsy, std::vector<double>& positionsz)
{
    tree->rescan_tree_v(positionsx, positionsy, positionsz, radii, volumes, gammas, ishydrogen);
}

//deposit current gammas on the overlap tree
void GaussVol::rescan_tree_gammas(void)
{
    tree->rescan_tree_g(gammas);
}

int GOverlap_Tree::nchildren_under_slot_r(int slot)
{
    int n = 0;
    if(overlaps[slot].children_count > 0) {
        n += overlaps[slot].children_count;
        //now calls itself on the children
        for(int i = 0; i < overlaps[slot].children_count; i++) {
            n += nchildren_under_slot_r(overlaps[slot].children_startindex + i);
        }
    }
    return n;
}

// returns number of overlaps for each atom 
void GaussVol::getstat(std::vector<int>& nov)
{
    nov.resize(natoms);
    for(int i = 0; i < natoms; i++) nov[i] = 0;
    for(int atom = 0; atom < natoms; atom++) {
        int slot = atom + 1;
        nov[atom] = tree->nchildren_under_slot_r(slot);
    }
}
