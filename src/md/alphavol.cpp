// Author: Moses KJ Chung
// Year:   2024

#include "alphavol.h"
#include "libfunc.h"
#include "mathConst.h"

namespace polmdqc
{
/////////////////////////////////////////
//                                     //
//  alphavol  --  Alpha volumes class  //
//                                     //
/////////////////////////////////////////

// "alphavol" computes the surface area, volume, mean, and
// gaussian curvature; also computes derivatives

template <bool compder>
void AlphaVol::alphavol(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::vector<Edge>& edges, std::vector<Face>& faces,
    real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
    real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord)
{
    int ia,ib,ic,id;
    int e1,e2,e3;
    int edge_list[6];
    real ra,ra2,rb,rb2,rc,rc2,rd,rd2;
    real rab,rac,rad,rbc,rbd,rcd;
    real rab2,rac2,rad2,rbc2,rbd2,rcd2;
    real val,val1S,val2S,val3S,val4S;
    real val1M1,val2M1,val3M1,val4M1;
    real val1G1;
    real d1,d2,d3,d4;
    real val1V,val2V,val3V,val4V;
    real coefval,coefvalS;
    real surfa,surfb,surfc;
    real gaussa,gaussb,gaussc;
    real vola,volb,volc,vold;
    real r,phi,l,dr,dphi,dl;
    real coefaS,coefbS,coefcS,coefdS;
    real coefaV,coefbV,coefcV,coefdV;
    real coefaM,coefbM,coefcM,coefdM;
    real coefaG,coefbG,coefcG,coefdG;
    real dsurfa2,dsurfb2;
    real dvola2,dvolb2;
    real u[3];
    real angle[6],cosine[6],sine[6];
    real deriv[6][6],deriv2[6][3],dg[3][3];
    real dsurfa3[3],dsurfb3[3],dsurfc3[3];
    real dvola3[3],dvolb3[3],dvolc3[3];
    real dvola[6],dvolb[6],dvolc[6],dvold[6];
    constexpr real twopi = 2 * pi;
    constexpr real coefe = 1.;
    constexpr real coeff = 2.;

    int nedges = edges.size();
    int nvertices = vertices.size();
    int nfaces = faces.size();
    int ntetra = tetra.size();

    // initialize results arrays
    for (int i = 0; i < nvertices; i++) {
        ballwsurf[i] = 0.;
        ballwvol[i] = 0.;
        ballwmean[i] = 0.;
        ballwgauss[i] = 0.;
    }

    // initialize edge and vertex info
    for (int i = 0; i < nedges; i++) {
        edges[i].gamma = 1.;
        edges[i].sigma = 1.;

        ia = edges[i].vertices[0];
        ib = edges[i].vertices[1];

        if (vertices[ia].status==0 || vertices[ib].status==0) continue;

        ra = vertices[ia].r; ra2 = ra*ra;
        rb = vertices[ib].r; rb2 = rb*rb;

        coefaS = vertices[ia].coefs; coefbS = vertices[ib].coefs;
        coefaV = vertices[ia].coefv; coefbV = vertices[ib].coefv;
        coefaM = vertices[ia].coefm; coefbM = vertices[ib].coefm;
        coefaG = vertices[ia].coefg; coefbG = vertices[ib].coefg;

        rab2 = dist2(vertices, ia, ib);
        rab = REAL_SQRT(rab2);

        twosph(ra, ra2, rb, rb2, rab, rab2, surfa, surfb, vola, volb, r, phi, l);

        edges[i].len = rab;
        edges[i].surf   = (coefaS*surfa + coefbS*surfb)/twopi;
        edges[i].vol    = (coefaV*vola + coefbV*volb)/twopi;
        edges[i].coefm1 = (coefaM*surfa/ra + coefbM*surfb/rb)/twopi;
        edges[i].coefm2 = coefe*(coefaM+coefbM)*r*phi/2;
        edges[i].coefg1 = (coefaG*surfa/ra2 + coefbG*surfb/rb2)/twopi;
        edges[i].coefg2 = coefe*(coefaG+coefbG)*l/2;
        edges[i].dsurf  = 0;
        edges[i].dvol   = 0;
        edges[i].dmean  = 0;
        edges[i].dgauss = 0;

    }

    for (int i = 0; i < nvertices; i++) {
        vertices[i].gamma = 1.;
    }

    // Contributions of four overlapping spheres:
    // We are using the weighted inclusion-exclusion formula:
    // Each tetrahedron in the Alpha Complex only contributes to the weight of each
    // its edges and each of its vertices

    for (int idx = 0; idx < ntetra; idx++) {
        if (tetra[idx].info[6]==0) continue;

        ia = tetra[idx].vertices[0];
        ib = tetra[idx].vertices[1];
        ic = tetra[idx].vertices[2];
        id = tetra[idx].vertices[3];

        if (vertices[ia].status==0 || vertices[ib].status==0
         || vertices[ic].status==0 || vertices[id].status==0) continue;

        ra = vertices[ia].r; ra2 = ra*ra;
        rb = vertices[ib].r; rb2 = rb*rb;
        rc = vertices[ic].r; rc2 = rc*rc;
        rd = vertices[id].r; rd2 = rd*rd;

        coefaS = vertices[ia].coefs; coefaV = vertices[ia].coefv;
        coefaM = vertices[ia].coefm; coefaG = vertices[ia].coefg;
        coefbS = vertices[ib].coefs; coefbV = vertices[ib].coefv;
        coefbM = vertices[ib].coefm; coefbG = vertices[ib].coefg;
        coefcS = vertices[ic].coefs; coefcV = vertices[ic].coefv;
        coefcM = vertices[ic].coefm; coefcG = vertices[ic].coefg;
        coefdS = vertices[id].coefs; coefdV = vertices[id].coefv;
        coefdM = vertices[id].coefm; coefdG = vertices[id].coefg;

        for (int iedge = 0; iedge < 6; iedge++) {
        // iedge is the edge number in the tetrahedron idx, with:
        // iedge = 1  (c,d)
        // iedge = 2  (b,d)
        // iedge = 3  (b,c)
        // iedge = 4  (a,d)
        // iedge = 5  (a,c)
        // iedge = 6  (a,b)
            edge_list[5-iedge] = tetra[idx].info_edge[iedge];
        }

        rab = edges[edge_list[0]].len; rab2 = rab*rab;
        rac = edges[edge_list[1]].len; rac2 = rac*rac;
        rad = edges[edge_list[2]].len; rad2 = rad*rad;
        rbc = edges[edge_list[3]].len; rbc2 = rbc*rbc;
        rbd = edges[edge_list[4]].len; rbd2 = rbd*rbd;
        rcd = edges[edge_list[5]].len; rcd2 = rcd*rcd;

        // characterize tetrahedron (A,B,C,D)
        tetdihedder<compder>(rab2, rac2, rad2, rbc2, rbd2, rcd2, angle, cosine, sine, deriv);

        // add fraction of tetrahedron that "belongs" to each ball
        tetvorder<compder>(ra2, rb2, rc2, rd2, rab, rac, rad, rbc,
        rbd, rcd, rab2, rac2, rad2, rbc2, rbd2, rcd2, cosine, sine,
        deriv, vola, volb, volc, vold, dvola, dvolb, dvolc, dvold);

        ballwvol[ia] += vola; ballwvol[ib] += volb;
        ballwvol[ic] += volc; ballwvol[id] += vold;

        if constexpr (compder) {
            for (int iedge = 0; iedge < 6; iedge++) {
                int i1 = edge_list[iedge];
                edges[i1].dvol += coefaV*dvola[iedge]
                +coefbV*dvolb[iedge]+coefcV*dvolc[iedge]
                +coefdV*dvold[iedge];
            }
        }

        // weights on each vertex: fraction of solid angle
        vertices[ia].gamma -= (angle[0]+angle[1]+angle[2])/2 - 0.250;
        vertices[ib].gamma -= (angle[0]+angle[3]+angle[4])/2 - 0.250;
        vertices[ic].gamma -= (angle[1]+angle[3]+angle[5])/2 - 0.250;
        vertices[id].gamma -= (angle[2]+angle[4]+angle[5])/2 - 0.250;

        // weights on each edge: fraction of dihedral angle
        for (int iedge = 0; iedge < 6; iedge++) {
            int i1 = edge_list[iedge];
            if (edges[i1].gamma != 0) {
                edges[i1].gamma -= angle[iedge];
                edges[i1].sigma -= angle[iedge];
            }
        }

        if constexpr (compder) {
            // Derivative: take into account the derivatives of the edge
            // weight in weightedinclusion-exclusion formula
            for (int iedge = 0; iedge < 6; iedge++) {
                int i1 = edge_list[iedge];
                val1S = edges[i1].surf;
                val1V = edges[i1].vol;
                val1M1 = edges[i1].coefm1 + edges[i1].coefm2;
                val1G1 = edges[i1].coefg1 + edges[i1].coefg2;
                for (int ie = 0; ie < 6; ie++) {
                    int je = edge_list[ie];
                    edges[je].dsurf  += val1S*deriv[iedge][ie];
                    edges[je].dvol   += val1V*deriv[iedge][ie];
                    edges[je].dmean  += val1M1*deriv[iedge][ie];
                    edges[je].dgauss += val1G1*deriv[iedge][ie];
                }
            }

            // Derivative: take into account the derivatives of the vertex
            // weight in weighted inclusion-exclusion formula

            val1S = ra2*coefaS; val2S = rb2*coefbS;
            val3S = rc2*coefcS; val4S = rd2*coefdS;

            val1V = ra2*ra*coefaV/3; val2V = rb2*rb*coefbV/3;
            val3V = rc2*rc*coefcV/3; val4V = rd2*rd*coefdV/3;

            val1M1 = ra*coefaM; val2M1 = rb*coefbM;
            val3M1 = rc*coefcM; val4M1 = rd*coefdM;

            for (int ie = 0; ie < 6; ie++) {
                int je = edge_list[ie];
                d1 = deriv[0][ie]+deriv[1][ie]+deriv[2][ie];
                d2 = deriv[0][ie]+deriv[3][ie]+deriv[4][ie];
                d3 = deriv[1][ie]+deriv[3][ie]+deriv[5][ie];
                d4 = deriv[2][ie]+deriv[4][ie]+deriv[5][ie];

                val = val1S*d1 + val2S*d2 + val3S*d3 + val4S*d4;
                edges[je].dsurf -= val;

                val = val1V*d1 + val2V*d2 + val3V*d3 + val4V*d4;
                edges[je].dvol -= val;

                val = val1M1*d1 + val2M1*d2 + val3M1*d3 + val4M1*d4;
                edges[je].dmean -= val;

                val = coefaG*d1 + coefbG*d2 + coefcG*d3 + coefdG*d4;
                edges[je].dgauss -= val;
            }
        }
    }

    // contribution of 3-balls (i.e. triangles of the alpha complex)

    for (int idx = 0; idx < nfaces; idx++) {
        coefval = faces[idx].gamma;
        if (coefval==0) continue;

        ia = faces[idx].vertices[0];
        ib = faces[idx].vertices[1];
        ic = faces[idx].vertices[2];

        if (vertices[ia].status==0 || vertices[ib].status==0
         || vertices[ic].status==0 ) continue;


        coefaS = vertices[ia].coefs; coefaV = vertices[ia].coefv;
        coefaM = vertices[ia].coefm; coefaG = vertices[ia].coefg;
        coefbS = vertices[ib].coefs; coefbV = vertices[ib].coefv;
        coefbM = vertices[ib].coefm; coefbG = vertices[ib].coefg;
        coefcS = vertices[ic].coefs; coefcV = vertices[ic].coefv;
        coefcM = vertices[ic].coefm; coefcG = vertices[ic].coefg;

        e1 = faces[idx].edges[0];
        e2 = faces[idx].edges[1];
        e3 = faces[idx].edges[2];

        ra = vertices[ia].r; ra2 = ra*ra;
        rb = vertices[ib].r; rb2 = rb*rb;
        rc = vertices[ic].r; rc2 = rc*rc;

        rab = edges[e1].len; rab2=rab*rab;
        rac = edges[e2].len; rac2=rac*rac;
        rbc = edges[e3].len; rbc2=rbc*rbc;

        threesphder<compder>(ra, rb, rc, ra2, rb2, rc2, rab, rac, rbc, rab2, rac2, rbc2,
        angle, deriv2, surfa, surfb, surfc, vola, volb, volc,
        dsurfa3, dsurfb3, dsurfc3, dvola3, dvolb3, dvolc3);

        threesphgss<compder>(ra, rb, rc, ra2, rb2, rc2, rab, rac, rbc,
        rab2, rac2, rbc2, gaussa, gaussb, gaussc, dg);

        ballwsurf[ia] += coefval*surfa;
        ballwsurf[ib] += coefval*surfb;
        ballwsurf[ic] += coefval*surfc;

        ballwvol[ia] += coefval*vola;
        ballwvol[ib] += coefval*volb;
        ballwvol[ic] += coefval*volc;

        ballwgauss[ia] += coeff*coefval*gaussa;
        ballwgauss[ib] += coeff*coefval*gaussb;
        ballwgauss[ic] += coeff*coefval*gaussc;

        edges[e1].sigma -= 2*coefval*angle[0];
        edges[e2].sigma -= 2*coefval*angle[1];
        edges[e3].sigma -= 2*coefval*angle[3];

        if constexpr (compder) {
            edges[e1].dsurf += coefval*(coefaS*dsurfa3[0]+coefbS*dsurfb3[0]+coefcS*dsurfc3[0]);
            edges[e2].dsurf += coefval*(coefaS*dsurfa3[1]+coefbS*dsurfb3[1]+coefcS*dsurfc3[1]);
            edges[e3].dsurf += coefval*(coefaS*dsurfa3[2]+coefbS*dsurfb3[2]+coefcS*dsurfc3[2]);

            edges[e1].dvol += coefval*(coefaV*dvola3[0]+coefbV*dvolb3[0]+coefcV*dvolc3[0]);
            edges[e2].dvol += coefval*(coefaV*dvola3[1]+coefbV*dvolb3[1]+coefcV*dvolc3[1]);
            edges[e3].dvol += coefval*(coefaV*dvola3[2]+coefbV*dvolb3[2]+coefcV*dvolc3[2]);

            edges[e1].dmean += coefval*(coefaM*dsurfa3[0]/ra+coefbM*dsurfb3[0]/rb+
                coefcM*dsurfc3[0]/rc);
            edges[e2].dmean += coefval*(coefaM*dsurfa3[1]/ra+coefbM*dsurfb3[1]/rb+
                coefcM*dsurfc3[1]/rc);
            edges[e3].dmean += coefval*(coefaM*dsurfa3[2]/ra+coefbM*dsurfb3[2]/rb+
                coefcM*dsurfc3[2]/rc);

            edges[e1].dmean += 2*coefval*(edges[e1].coefm2*deriv2[0][0]+
                edges[e2].coefm2*deriv2[1][0]+edges[e3].coefm2*deriv2[3][0]);
            edges[e2].dmean += 2*coefval*(edges[e1].coefm2*deriv2[0][1]+
                edges[e2].coefm2*deriv2[1][1]+edges[e3].coefm2*deriv2[3][1]);
            edges[e3].dmean += 2*coefval*(edges[e1].coefm2*deriv2[0][2]+
                edges[e2].coefm2*deriv2[1][2]+edges[e3].coefm2*deriv2[3][2]);

            edges[e1].dgauss += coefval*(coefaG*dsurfa3[0]/ra2+coefbG*dsurfb3[0]/rb2+
                coefcG*dsurfc3[0]/rc2);
            edges[e2].dgauss += coefval*(coefaG*dsurfa3[1]/ra2+coefbG*dsurfb3[1]/rb2+
                coefcG*dsurfc3[1]/rc2);
            edges[e3].dgauss += coefval*(coefaG*dsurfa3[2]/ra2+coefbG*dsurfb3[2]/rb2+
                coefcG*dsurfc3[2]/rc2);

            edges[e1].dgauss += 2*coefval*(edges[e1].coefg2*deriv2[0][0]+
                edges[e2].coefg2*deriv2[1][0]+edges[e3].coefg2*deriv2[3][0]);
            edges[e2].dgauss += 2*coefval*(edges[e1].coefg2*deriv2[0][1]+
                edges[e2].coefg2*deriv2[1][1]+edges[e3].coefg2*deriv2[3][1]);
            edges[e3].dgauss += 2*coefval*(edges[e1].coefg2*deriv2[0][2]+
                edges[e2].coefg2*deriv2[1][2]+edges[e3].coefg2*deriv2[3][2]);

            edges[e1].dgauss += coeff*coefval*(coefaG*dg[0][0]+coefbG*dg[1][0]+coefcG*dg[2][0]);
            edges[e2].dgauss += coeff*coefval*(coefaG*dg[0][1]+coefbG*dg[1][1]+coefcG*dg[2][1]);
            edges[e3].dgauss += coeff*coefval*(coefaG*dg[0][2]+coefbG*dg[1][2]+coefcG*dg[2][2]);
        }
    }

    // now add contribution of two-sphere
    real coefeps = 1.e-10;
    for (int iedge = 0; iedge < nedges; iedge++) {
        coefval = edges[iedge].gamma;
        coefvalS = edges[iedge].sigma;

        if (REAL_ABS(coefval) < coefeps) continue;

        ia = edges[iedge].vertices[0];
        ib = edges[iedge].vertices[1];

        if (vertices[ia].status==0 || vertices[ib].status==0) continue;

        coefaS = vertices[ia].coefs; coefaV = vertices[ia].coefv;
        coefaM = vertices[ia].coefm; coefaG = vertices[ia].coefg;
        coefbS = vertices[ib].coefs; coefbV = vertices[ib].coefv;
        coefbM = vertices[ib].coefm; coefbG = vertices[ib].coefg;

        ra = vertices[ia].r; ra2 = ra*ra;
        rb = vertices[ib].r; rb2 = rb*rb;

        rab = edges[iedge].len; rab2 = rab*rab;

        twosphder<compder>(ra, ra2, rb, rb2, rab, rab2, surfa, surfb,
        vola, volb, r, phi, l, dsurfa2, dsurfb2, dvola2, dvolb2, dr, dphi, dl);

        ballwsurf[ia] -= coefval*surfa;
        ballwsurf[ib] -= coefval*surfb;
        ballwvol[ia]  -= coefval*vola;
        ballwvol[ib]  -= coefval*volb;

        val = coefe*pi*coefvalS*r*phi;
        ballwmean[ia] -= val;
        ballwmean[ib] -= val;

        val = coefe*pi*coefvalS*l;
        ballwgauss[ia] -= val;
        ballwgauss[ib] -= val;

        if constexpr (compder) {
            edges[iedge].dsurf  -= coefval* (coefaS*dsurfa2 + coefbS*dsurfb2);
            edges[iedge].dvol   -= coefval* (coefaV*dvola2 + coefbV*dvolb2);
            edges[iedge].dmean  -= coefval* (coefaM*dsurfa2/ra + coefbM*dsurfb2/rb);
            edges[iedge].dmean  -= coefe*coefvalS*pi*(r*dphi+phi*dr)*(coefaM+coefbM);
            edges[iedge].dgauss -= coefval* (coefaG*dsurfa2/ra2 + coefbG*dsurfb2/rb2);
            edges[iedge].dgauss -= coefe*coefvalS*pi*dl*(coefaG+coefbG);
        }
    }

    // now loop over vertices
    for (int i = 4; i < nvertices; i++) {
        coefval = vertices[i].gamma;
        if (vertices[i].info[0]==0) continue;
        if (vertices[i].info[7]==0) continue;
        if (coefval==0) continue;
        if (vertices[i].status==0) continue;

        ra = vertices[i].r; ra2 = ra*ra;
        surfa = 4*pi*ra*ra;
        vola  = surfa*ra/3;
        ballwsurf[i]  += coefval*surfa;
        ballwvol[i]   += coefval*vola;
        if (ra2 > 0) {
            ballwmean[i]  += ballwsurf[i]/ra;
            ballwgauss[i] += ballwsurf[i]/ra2;
        }
    }

    // compute total surface, volume (weighted, and unweighted)
    for (int i = 4; i < nvertices; i++) {
        if (vertices[i].info[0]==0) continue;
        if (vertices[i].status==0) continue;

        coefaS = vertices[i].coefs; coefaV = vertices[i].coefv;
        coefaM = vertices[i].coefm; coefaG = vertices[i].coefg;

        ballwsurf[i] *= coefaS;
        ballwvol[i] *= coefaV;
        ballwmean[i] *= coefaM;
        ballwgauss[i] *= coefaG;
    }

    // shift as 4 first vertices are pseudo atoms
    int nballs = 0;
    for (int i = 0; i < nvertices; i++) if (vertices[i].status==1) nballs++;
    for (int i = 0; i < nballs; i++) {
        ballwsurf[i]   = ballwsurf[i+4];
        ballwvol[i]    = ballwvol[i+4];
        ballwmean[i]   = ballwmean[i+4];
        ballwgauss[i]  = ballwgauss[i+4];
    }

    if constexpr (!compder) return;

    // convert derivatives wrt to distance to derivatives wrt to coordinates
    for (int i = 0; i < 3*nvertices; i++) {
        dsurf_coord[i] = 0.;
        dvol_coord[i] = 0.;
        dmean_coord[i] = 0.;
        dgauss_coord[i] = 0.;
    }

    for (int iedge = 0; iedge < nedges; iedge++) {
        ia = edges[iedge].vertices[0];
        ib = edges[iedge].vertices[1];

        for (int i = 0; i < 3; i++) {
            u[i] = vertices[ia].coord[i] - vertices[ib].coord[i];
        }

        rab  = edges[iedge].len;
        val1S  = edges[iedge].dsurf/rab;
        val1V  = edges[iedge].dvol/rab;
        val1M1 = edges[iedge].dmean/rab;
        val1G1 = edges[iedge].dgauss/rab;

        for (int j = 0; j < 3; j++) {
            dsurf_coord[3*ia+j]  += u[j]*val1S;
            dsurf_coord[3*ib+j]  -= u[j]*val1S;
            dvol_coord[3*ia+j]   += u[j]*val1V;
            dvol_coord[3*ib+j]   -= u[j]*val1V;
            dmean_coord[3*ia+j]  += u[j]*val1M1;
            dmean_coord[3*ib+j]  -= u[j]*val1M1;
            dgauss_coord[3*ia+j] += u[j]*val1G1;
            dgauss_coord[3*ib+j] -= u[j]*val1G1;
        }
    }

    // shift as 4 first vertices are pseudo atoms
    for (int i = 0; i < 3*nballs; i++) {
        dsurf_coord[i]   = dsurf_coord[i+12];
        dvol_coord[i]    = dvol_coord[i+12];
        dmean_coord[i]   = dmean_coord[i+12];
        dgauss_coord[i]  = dgauss_coord[i+12];
    }
}

inline real AlphaVol::dist2(std::vector<Vertex>& vertices, int n1, int n2)
{
    real x;
    real dist = 0;
    for(int i = 0; i < 3; i++) {
        x = vertices[n1].coord[i] - vertices[n2].coord[i];
        dist += x*x;
    }

    return dist;
}

/////////////////////////////////////////////////////////
//                                                     //
//  twosph  --  surf & vol of two overlapping spheres  //
//                                                     //
/////////////////////////////////////////////////////////

// "twosph" calculates the surface area and volume of the
// intersection of two spheres

inline void AlphaVol::twosph(real ra, real ra2, real rb, real rb2,
    real rab, real rab2, real& surfa, real& surfb,
    real& vola, real& volb, real& r, real& phi, real& l)
{
    real cosine,vala,valb,lambda,ha,hb;
    real Aab,sa,ca,sb,cb;
    constexpr real twopi = 2 * pi;

    // Get distance between center of sphere A and Voronoi plane
    // between A and B
    lambda = plane_dist(ra2, rb2, rab2);
    valb = lambda*rab;
    vala = rab-valb;

    // Get height of the cap of sphere A occluded by sphere B
    ha = ra - vala;

    // same for sphere B ...
    hb = rb - valb;

    // get surfaces of intersection
    surfa = twopi*ra*ha;
    surfb = twopi*rb*hb;

    // now get volume
    Aab = pi*(ra2-vala*vala);

    sa = ra*(surfa);
    ca = vala*Aab;

    vola = (sa-ca)/3;

    sb = rb*(surfb);
    cb = valb*Aab;

    volb = (sb-cb)/3;

    // get radius of the circle of intersection between the two spheres
    r = REAL_SQRT(ra2 - vala*vala);

    // get angle between normals of the sphere at a point on this circle
    cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
    if (REAL_ABS(cosine - 1) < eps) cosine = 1;
    else if (REAL_ABS(cosine + 1) < eps) cosine = -1;
    phi = REAL_ACOS(cosine);

    l = vala/ra + valb/rb;
}

///////////////////////////////////////////////////////////////////
//                                                               //
//  twosphder  --  surf & vol derivs of two overlapping spheres  //
//                                                               //
///////////////////////////////////////////////////////////////////

// "twosphder" calculates the surface area and volume derivatives
// of the intersection of two spheres

template <bool compder>
inline void AlphaVol::twosphder(real ra, real ra2, real rb, real rb2, real rab, real rab2,
    real& surfa, real& surfb, real& vola, real& volb, real& r, real& phi, real& l,
    real& dsurfa, real& dsurfb, real& dvola, real& dvolb, real& dr, real& dphi, real& dl)
{
    real cosine,vala,valb,lambda,ha,hb;
    real Aab,sa,ca,sb,cb;
    real dera,derb;
    constexpr real twopi = 2 * pi;

    // Get distance between center of sphere A and Voronoi plane
    // between A and B
    lambda = plane_dist(ra2, rb2, rab2);
    valb = lambda*rab;
    vala = rab-valb;

    // get height of the cap of sphere A occluded by sphere B
    ha = ra - vala;

    // same for sphere B ...
    hb = rb - valb;

    // get surfaces of intersection
    surfa = twopi*ra*ha;
    surfb = twopi*rb*hb;

    // now get volume
    Aab = pi*(ra2-vala*vala);

    sa = ra*(surfa);
    ca = vala*Aab;

    vola = (sa-ca)/3;

    sb = rb*(surfb);
    cb = valb*Aab;

    volb = (sb-cb)/3;

    // get radius of the circle of intersection between the two spheres
    r = REAL_SQRT(ra2 - vala*vala);

    // get angle between normals of the sphere at a point on this circle
    cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
    if (REAL_ABS(cosine - 1) < eps) cosine = 1;
    else if (REAL_ABS(cosine + 1) < eps) cosine = -1;
    phi = REAL_ACOS(cosine);
    l = vala/ra + valb/rb;

    if constexpr (!compder) return;

    dera = - lambda;
    derb = lambda - 1;

    dsurfa = twopi*ra*dera;
    dsurfb = twopi*rb*derb;

    dvola = -Aab*lambda;
    dvolb = -(dvola) - Aab;

    dr   = -vala*lambda/(r);
    dphi = rab/(ra*rb*REAL_SQRT(1-cosine*cosine));
    dl   = lambda/ra + (1.0-lambda)/rb;
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  threesphder  --  surf & vol derivs of three overlapping spheres  //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// "threesphder" calculates the surface area and volume derivatives
// of the intersection of three spheres

template <bool compder>
inline void AlphaVol::threesphder(real ra, real rb,real rc, real ra2,
    real rb2, real rc2, real rab, real rac, real rbc,
    real rab2, real rac2, real rbc2, real *angle, real deriv[6][3],
    real& surfa, real& surfb, real& surfc, real& vola, real& volb, real& volc,
    real* dsurfa, real* dsurfb, real* dsurfc, real* dvola, real* dvolb, real* dvolc)
{
    real a1,a2,a3,s2,c1,c2;
    real seg_ang_ab,seg_ang_ac,seg_ang_bc;
    real ang_dih_ap,ang_dih_bp,ang_dih_cp;
    real val1,val2,val3,l1,l2,l3;
    real val1b,val2b,val3b;
    real ang_abc,ang_acb,ang_bca;
    real cos_abc,cos_acb,cos_bca;
    real sin_abc,sin_acb,sin_bca;
    real s_abc,s_acb,s_bca;
    real rho_ab2,rho_ac2,rho_bc2;
    real drho_ab2,drho_ac2,drho_bc2;
    real val_abc,val_acb,val_bca;
    real val2_abc,val2_acb,val2_bca;
    real der_val1b,der_val1,der_val2b,der_val2,der_val3b,der_val3;
    real cosine[6],sine[6];
    constexpr real twopi = 2 * pi;

    l1 = plane_dist(ra2, rb2, rab2);
    l2 = plane_dist(ra2, rc2, rac2);
    l3 = plane_dist(rb2, rc2, rbc2);

    val1 = l1*rab; val2 = l2*rac; val3 = l3*rbc;
    val1b = rab - val1; val2b = rac - val2; val3b = rbc - val3;

    // We consider the tetrahedron (A,B,C,P) where P is the
    // point of intersection of the three spheres such that (A,B,C,P) is ccw.
    // The edge lengths in this tetrahedron are: rab, rac, rAP=ra, rbc, rBP=rb, rCP=rc

    tetdihedder3<compder>(rab2, rac2, ra2, rbc2, rb2, rc2, angle, cosine, sine, deriv);

    // the seg_ang_ are the dihedral angles around the three edges AB, AC and BC

    seg_ang_ab = angle[0];
    seg_ang_ac = angle[1];
    seg_ang_bc = angle[3];

    // the ang_dih_ are the dihedral angles around the three edges AP, BP and CP
    ang_dih_ap = angle[2];
    ang_dih_bp = angle[4];
    ang_dih_cp = angle[5];

    a1 = ra*(1-2*ang_dih_ap);
    a2 = 2*seg_ang_ab*val1b;
    a3 = 2*seg_ang_ac*val2b;

    surfa = twopi*ra*(a1 - a2 - a3);

    a1 = rb*(1-2*ang_dih_bp);
    a2 = 2*seg_ang_ab*val1;
    a3 = 2*seg_ang_bc*val3b;

    surfb = twopi*rb*(a1 - a2 - a3);

    a1 = rc*(1-2*ang_dih_cp);
    a2 = 2*seg_ang_ac*val2;
    a3 = 2*seg_ang_bc*val3;

    surfc = twopi*rc*(a1 - a2 - a3);

    // compute volumes of the three caps
    ang_abc = twopi*seg_ang_ab;
    ang_acb = twopi*seg_ang_ac;
    ang_bca = twopi*seg_ang_bc;

    cos_abc = cosine[0];
    sin_abc = sine[0];
    cos_acb = cosine[1];
    sin_acb = sine[1];
    cos_bca = cosine[3];
    sin_bca = sine[3];

    rho_ab2 = ra2 - val1b*val1b;
    rho_ac2 = ra2 - val2b*val2b;
    rho_bc2 = rb2 - val3b*val3b;

    val_abc = ang_abc - sin_abc*cos_abc; s_abc = rho_ab2*val_abc;
    val_acb = ang_acb - sin_acb*cos_acb; s_acb = rho_ac2*val_acb;
    val_bca = ang_bca - sin_bca*cos_bca; s_bca = rho_bc2*val_bca;

    s2 = ra*(surfa);
    c1 = val1b*s_abc;
    c2 = val2b*s_acb;

    vola = (s2 - c1 - c2)/3;

    s2 = rb*(surfb);
    c1 = val1*s_abc;
    c2 = val3b*s_bca;

    volb = (s2 - c1 - c2)/3;

    s2 = rc*(surfc);
    c1 = val2*s_acb;
    c2 = val3*s_bca;

    volc = (s2 - c1 - c2)/3;

    if constexpr (!compder) return;

    der_val1b = l1; der_val1  = 1-l1;
    der_val2b = l2; der_val2  = 1-l2;
    der_val3b = l3; der_val3  = 1-l3;

    dsurfa[0] = -2*ra*(
        twopi*seg_ang_ab*der_val1b +
        (ra*deriv[2][0] +
        val1b*deriv[0][0] +val2b*deriv[1][0]));
    dsurfa[1] = -2*ra*(
        twopi*seg_ang_ac*der_val2b +
        (ra*deriv[2][1] +
        val1b*deriv[0][1] +val2b*deriv[1][1]));
    dsurfa[2] = -2*ra*( ra*deriv[2][2] +
        val1b*deriv[0][2]+val2b*deriv[1][2]);

    dsurfb[0] = -2*rb*(
        twopi*seg_ang_ab*der_val1
        +(rb*deriv[4][0]+
        val1*deriv[0][0]+val3b*deriv[3][0]));
    dsurfb[1] = -2*rb*(rb*deriv[4][1]+
        val1*deriv[0][1]+val3b*deriv[3][1]);
    dsurfb[2] = -2*rb*(
        twopi*seg_ang_bc*der_val3b
        +(rb*deriv[4][2]+
        val1*deriv[0][2]+val3b*deriv[3][2]));

    dsurfc[0] = -2*rc*(rc*deriv[5][0]+
            val2*deriv[1][0]+val3*deriv[3][0]);
    dsurfc[1] = -2*rc*(
        twopi*seg_ang_ac*der_val2
        +(rc*deriv[5][1]+
        val2*deriv[1][1]+val3*deriv[3][1]));
    dsurfc[2] = -2*rc*(
        twopi*seg_ang_bc*der_val3
        +(rc*deriv[5][2]+
        val2*deriv[1][2]+val3*deriv[3][2]));

    drho_ab2 = -2*der_val1b*val1b;
    drho_ac2 = -2*der_val2b*val2b;
    drho_bc2 = -2*der_val3b*val3b;

    val2_abc = rho_ab2*(1 - cos_abc*cos_abc + sin_abc*sin_abc);
    val2_acb = rho_ac2*(1 - cos_acb*cos_acb + sin_acb*sin_acb);
    val2_bca = rho_bc2*(1 - cos_bca*cos_bca + sin_bca*sin_bca);

    dvola[0] = ra*dsurfa[0] - der_val1b*s_abc -
        (val1b*deriv[0][0]*val2_abc + val2b*deriv[1][0]*val2_acb)
        - val1b*drho_ab2*val_abc;
    dvola[0] = dvola[0]/3;
    dvola[1] = ra*dsurfa[1] - der_val2b*s_acb -
        (val1b*deriv[0][1]*val2_abc + val2b*deriv[1][1]*val2_acb)
        - val2b*drho_ac2*val_acb;
    dvola[1] = dvola[1]/3;
    dvola[2] = ra*dsurfa[2] -
        (val1b*deriv[0][2]*val2_abc + val2b*deriv[1][2]*val2_acb);
    dvola[2] = dvola[2]/3;

    dvolb[0] = rb*dsurfb[0] - der_val1*s_abc -
        (val1*deriv[0][0]*val2_abc + val3b*deriv[3][0]*val2_bca)
        - val1*drho_ab2*val_abc;
    dvolb[0] = dvolb[0]/3;
    dvolb[1] = rb*dsurfb[1] -
        (val1*deriv[0][1]*val2_abc + val3b*deriv[3][1]*val2_bca);
    dvolb[1] = dvolb[1]/3;
    dvolb[2] = rb*dsurfb[2] - der_val3b*s_bca -
        (val1*deriv[0][2]*val2_abc + val3b*deriv[3][2]*val2_bca)
        - val3b*drho_bc2*val_bca;
    dvolb[2] = dvolb[2]/3;

    dvolc[0] = rc*dsurfc[0] -
        (val2*deriv[1][0]*val2_acb + val3*deriv[3][0]*val2_bca);
    dvolc[0] = dvolc[0]/3;
    dvolc[1] = rc*dsurfc[1] - der_val2*s_acb -
        (val2*deriv[1][1]*val2_acb + val3*deriv[3][1]*val2_bca)
        - val2*drho_ac2*val_acb;
    dvolc[1] = dvolc[1]/3;
    dvolc[2] = rc*dsurfc[2] - der_val3*s_bca -
        (val2*deriv[1][2]*val2_acb + val3*deriv[3][2]*val2_bca)
        - val3*drho_bc2*val_bca;
    dvolc[2] = dvolc[2]/3;
}

inline real AlphaVol::plane_dist(real ra2, real rb2, real rab2)
{
    real lambda = 0.50 - (ra2-rb2)/(2*rab2);

    return lambda;
}

////////////////////////////////////////////////////////////////
//                                                            //
//  tetdihedder  --  tetrahedron dihedral angles derivatives  //
//                                                            //
////////////////////////////////////////////////////////////////

// "tetdihedder" computes the derivative of the six
// dihedral angles of atetrahedronfrom its edge lengths

template <bool compder>
inline void AlphaVol::tetdihedder(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq, real r34sq, real* angle,
    real* cosine, real* sine, real deriv[6][6])
{
    real val1,val2,val3,val4,vala;
    real val123,val124,val134,val234;
    real val213,val214,val314,val324,val312;
    real det12,det13,det14,det23,det24,det34;
    real minori[4];
    real dminori[4][6] = {0};
    real det[6],dnum[6][6],val[4];
    real dist[6];
    constexpr real twopi = 2 * pi;

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    if constexpr (compder) val[0] = val4; val[1] = val3; val[2] = val2; val[3] = val1;

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij)
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;
    det14 = -2*r14sq*val123 - val124*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;
    val312 = r12sq -r13sq -r23sq;
    val314 = r14sq -r13sq -r34sq;
    val324 = r24sq -r23sq -r34sq;

    det23 = -2*r23sq*val214 - val213*val234;
    det24 = -2*r24sq*val213 - val214*val234;
    det34 = -2*r34sq*val312 - val314*val324;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det14*val2*val3;
    cosine[3] = det23*val1*val4;
    cosine[4] = det24*val2*val4;
    cosine[5] = det34*val3*val4;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(cosine[i] - 1) < eps) cosine[i] = 1;
        else if (REAL_ABS(cosine[i] + 1) < eps) cosine[i] = -1;
    }

    for (int i = 0; i < 6; i++) {
        angle[i] = REAL_ACOS(cosine[i]);
        sine[i]  = REAL_SIN(angle[i]);
        angle[i] /= twopi;
    }

    if constexpr (!compder) return;

    det[5] = det12; det[4] = det13; det[3] = det14;
    det[2] = det23; det[1] = det24; det[0] = det34;

    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r14sq);
    dist[3] = REAL_SQRT(r23sq); dist[4] = REAL_SQRT(r24sq); dist[5] = REAL_SQRT(r34sq);

    // Now compute derivatives of the angles with respect to the edge lengths
    // Since (see above):
    //                       num(i,j)
    // cos(ang(i,j)) = --------------------
    //                  sqrt(M(i,i)*M(j,j))

    //   d(ang(i,j))                       dnum(i,j)                             (M(i,i)dM(j,j) +M(j,j)*dM(i,i))
    // ------------sin(ang(i,j)) =  -----------------------   -0.5*num(i,j) -----------------------------------------
    //    dr(a,b)                   sqrt(M(i,i)M(j,j))dr(a,b)                    M(i,i)M(j,j) sqrt(M(i,i)M(j,j))

    // which we can rewrite as:

    //   d(ang(i,j))                cosine(i,j)   dnum(i,j)                    dM(j,j) +  dM(i,i))
    // ------------sin(ang(i,j)) = -----------  -----------  -0.5*cosine(i,j)( -------- + ---------)
    //    dr(a,b)                    num(i,j)     dr(a,b)                       M(j,j)      M(i,i)

    dminori[0][3] = -(val234 + 2*r24sq); dminori[0][4] = -(val234 + 2*r23sq); dminori[0][5] = val234;
    dminori[1][1] = -(val134 + 2*r14sq); dminori[1][2] = -(val134 + 2*r13sq); dminori[1][5] = val134;
    dminori[2][0] = -(val124 + 2*r14sq); dminori[2][2] = -(val124 + 2*r12sq); dminori[2][4] = val124;
    dminori[3][0] = -(val123 + 2*r13sq); dminori[3][1] = -(val123 + 2*r12sq); dminori[3][3] = val123;

    dnum[5][0] = -2*val134+val123+val124; dnum[5][1] = 2*r12sq + val124; dnum[5][2] = 2*r12sq + val123;
    dnum[5][3] = -val124; dnum[5][4] = -val123; dnum[5][5] = -2*r12sq;

    dnum[4][0] = 2*r13sq+val134; dnum[4][1] = -2*val124 + val123 + val134; dnum[4][2] = 2*r13sq + val123;
    dnum[4][3] = -val134; dnum[4][4] = -2*r13sq; dnum[4][5] = -val123;

    dnum[3][0] = 2*r14sq+val134; dnum[3][1] = 2*r14sq + val124; dnum[3][2] = -2*val123 + val124 + val134;
    dnum[3][3] = -2*r14sq; dnum[3][4] = -val134; dnum[3][5] = -val124;

    dnum[2][0] = 2*r23sq+val234; dnum[2][1] = -val234; dnum[2][2] = -2*r23sq;
    dnum[2][3] = -2*val214+val213+val234; dnum[2][4] = 2*r23sq+val213; dnum[2][5] = -val213;

    dnum[1][0] = 2*r24sq+val234; dnum[1][1] = -2*r24sq; dnum[1][2] = -val234;
    dnum[1][3] = 2*r24sq+val214; dnum[1][4] = -2*val213 + val214 + val234; dnum[1][5] = -val214;

    dnum[0][0] = -2*r34sq; dnum[0][1] = 2*r34sq+val324; dnum[0][2] = -val324;
    dnum[0][3] = 2*r34sq+val314; dnum[0][4] = -val314; dnum[0][5] = -2*val312 + val314 + val324;

    int k = 0;
    int jj;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            jj = 5-k;
            if (det[k] != 0) {
                vala = cosine[jj]/sine[jj];
                val1 = -vala/det[k];
                val2 = vala/minori[j];
                val3 = vala/minori[i];
                for (int l = 0; l < 6; l++) {
                    deriv[jj][l] = val1*dnum[k][l]+val2*dminori[j][l]+val3*dminori[i][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            else {
                vala = -val[i]*val[j]/sine[jj];
                for (int l = 0; l < 6; l++) {
                    deriv[jj][l] = vala*dnum[k][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            k++;
        }
    }
}

/////////////////////////////////////////////////////////////////
//                                                             //
//  tetdihedder3  --  compute dihedral angles and derivatives  //
//                                                             //
/////////////////////////////////////////////////////////////////

// "tetdihedder3" computes the six dihedral angles of a tetrahedron
// ABCD from its edge lengths as well as their derivatives with
// respect to the 3 edge lengths AB, AC and BC

template <bool compder>
inline void AlphaVol::tetdihedder3(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq, real r34sq, real* angle,
    real* cosine, real* sine, real deriv[6][3])
{
    real val1,val2,val3,val4,vala;
    real val123,val124,val134,val234;
    real val213,val214,val314,val324,val312;
    real det12,det13,det14,det23,det24,det34;
    real minori[4];
    real dminori[4][3] = {0};
    real dist[3],det[6],dnum[6][3],val[4];
    constexpr real twopi = 2 * pi;

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    val[0] = val4; val[1] = val3; val[2] = val2; val[3] = val1;

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij)
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;
    det14 = -2*r14sq*val123 - val124*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;
    val312 = r12sq -r13sq -r23sq;
    val314 = r14sq -r13sq -r34sq;
    val324 = r24sq -r23sq -r34sq;

    det23 = -2*r23sq*val214 - val213*val234;
    det24 = -2*r24sq*val213 - val214*val234;
    det34 = -2*r34sq*val312 - val314*val324;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det14*val2*val3;
    cosine[3] = det23*val1*val4;
    cosine[4] = det24*val2*val4;
    cosine[5] = det34*val3*val4;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(cosine[i] - 1) < eps) cosine[i] = 1;
        else if (REAL_ABS(cosine[i] + 1) < eps) cosine[i] = -1;
    }

    for (int i = 0; i < 6; i++) {
        angle[i] = REAL_ACOS(cosine[i]);
        sine[i]  = REAL_SIN(angle[i]);
        angle[i] /= twopi;
    }

    if constexpr (!compder) return;

    // Now compute derivatives of the angles with respect to the edge lengths
    // Since (see above):
    //                      num(i,j)
    // cos(ang(i,j)) = --------------------
    //                 sqrt(M(i,i)*M(j,j))

    //  d(ang(i,j))                      dnum(i,j)                             (M(i,i)dM(j,j) +M(j,j)*dM(i,i))
    // ------------sin(ang(i,j)) =  -----------------------   -0.5*num(i,j) -----------------------------------------
    //    dr(a,b)                   sqrt(M(i,i)M(j,j))dr(a,b)                  M(i,i)M(j,j) sqrt(M(i,i)M(j,j))

    // which we can rewrite as:

    // d(ang(i,j))                 cosine(i,j)    dnum(i,j)                    dM(j,j) +  dM(i,i))
    // ------------sin(ang(i,j)) = -----------  -----------  -0.5*cosine(i,j)( -------- + ---------)
    //     dr(a,b)                   num(i,j)       dr(a,b)                      M(j,j)      M(i,i)

    det[5] = det12; det[4] = det13; det[3] = det14;
    det[2] = det23; det[1] = det24; det[0] = det34;
    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r23sq);

    dminori[0][2] = -(val234 + 2*r24sq);
    dminori[1][1] = -(val134 + 2*r14sq);
    dminori[2][0] = -(val124 + 2*r14sq);
    dminori[3][0] = -(val123 + 2*r13sq); dminori[3][1] = -(val123 + 2*r12sq); dminori[3][2] = val123;

    dnum[5][0] = -2*val134+val123+val124; dnum[5][1] = 2*r12sq + val124; dnum[5][2] = -val124;
    dnum[4][0] = 2*r13sq+val134; dnum[4][1] = -2*val124 + val123 + val134; dnum[4][2] = -val134;
    dnum[3][0] = 2*r14sq+val134; dnum[3][1] = 2*r14sq + val124; dnum[3][2] = -2*r14sq;
    dnum[2][0] = 2*r23sq+val234; dnum[2][1] = -val234; dnum[2][2] = -2*val214 + val213 + val234;
    dnum[1][0] = 2*r24sq+val234; dnum[1][1] = -2*r24sq; dnum[1][2] = 2*r24sq + val214;
    dnum[0][0] = -2*r34sq; dnum[0][1] = 2*r34sq+val324; dnum[0][2] = 2*r34sq + val314;

    int k = 0;
    int jj;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            jj = 5-k;
            if (det[k] != 0) {
                vala = cosine[jj]/sine[jj];
                val1 = -vala/det[k];
                val2 = vala/minori[j];
                val3 = vala/minori[i];
                for (int l = 0; l < 3; l++) {
                    deriv[jj][l] = val1*dnum[k][l]
                    +val2*dminori[j][l]+val3*dminori[i][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            else {
                vala = -val[i]*val[j]/sine[jj];
                for (int l = 0; l < 3; l++) {
                    deriv[jj][l] = vala*dnum[k][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            k++;
        }
    }
}

///////////////////////////////////////////////////////
//                                                   //
//  tet3dihedcos  --  compute three of siz dihedral  //
//                                                   //
///////////////////////////////////////////////////////

// "tet3dihedcos" computes three of the six dihedral angles
// of a tetrahedron from its edge lengths

template <bool compder>
inline void AlphaVol::tet3dihedcos(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq,real r34sq, real* cosine, real deriv[3][3])
{
    real val1, val2, val3, val4;
    real val123, val124, val134, val234;
    real val213, val214;
    real det12, det13, det23;
    real minori[4];
    real dminori[4][3] = {0};
    real dnum[3][3];
    real dist[3];

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij)
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;

    det23 = -2*r23sq*val214 - val213*val234;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det23*val1*val4;

    if constexpr (!compder) return;

    dminori[0][2] = -(val234 + 2*r24sq);
    dminori[1][1] = -(val134 + 2*r14sq);
    dminori[2][0] = -(val124 + 2*r14sq);
    dminori[3][0] = -(val123 + 2*r13sq);
    dminori[3][1] = -(val123 + 2*r12sq);
    dminori[3][2] = val123;

    dnum[0][0] = -2*val134+val123+val124;
    dnum[0][1] = 2*r12sq + val124;
    dnum[0][2] = -val124;

    dnum[1][0] = 2*r13sq + val134;
    dnum[1][1] = -2*val124 + val123 + val134;
    dnum[1][2] = -val134;

    dnum[2][0] = 2*r23sq + val234;
    dnum[2][1] = -val234;
    dnum[2][2] = -2*val214 + val213 + val234;

    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r23sq);

    for (int i = 0; i < 3; i++) {
        deriv[0][i] = dnum[0][i]*val1*val2 - cosine[0]*
        (dminori[2][i]/minori[2] + dminori[3][i]/minori[3]);
        deriv[1][i] = dnum[1][i]*val1*val3 - cosine[1]*
        (dminori[1][i]/minori[1] + dminori[3][i]/minori[3]);
        deriv[2][i] = dnum[2][i]*val1*val4 - cosine[2]*
        (dminori[0][i]/minori[0] + dminori[3][i]/minori[3]);
        deriv[0][i] *= 2*dist[i];
        deriv[1][i] *= 2*dist[i];
        deriv[2][i] *= 2*dist[i];
    }
}

/////////////////////////////////////////////////////////
//                                                     //
//  tetvorder  --  compute intersection of four balls  //
//                                                     //
/////////////////////////////////////////////////////////

// "tetvorder" computes the volume of the intersection of
// the tetrahedron formed by the center of 4 balls with
// the Voronoi cells corresponding to these balls

template <bool compder>
inline void AlphaVol::tetvorder(real ra2,real rb2,real rc2,real rd2,
    real rab, real rac, real rad, real rbc, real rbd,
    real rcd, real rab2, real rac2, real rad2,real rbc2,
    real rbd2, real rcd2, real* cos_ang, real* sin_ang,
    real deriv[6][6], real& vola, real& volb, real& volc,
    real& vold, real* dvola, real* dvolb, real* dvolc, real* dvold)
{
    real l1,l2,l3,l4,l5,l6;
    real val1,val2,val3,val4,val5,val6;
    real val1b,val2b,val3b,val4b,val5b,val6b;
    real val_ab,val_ac,val_bc,val_ad,val_bd,val_cd;
    real val1_ab,val1_ac,val1_ad,val1_bc,val1_bd,val1_cd;
    real val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd;
    real cos_abc,cos_acb,cos_bca,cos_abd,cos_adb,cos_bda;
    real cos_acd,cos_adc,cos_cda,cos_bcd,cos_bdc,cos_cdb;
    real rho_ab2,rho_ac2,rho_ad2,rho_bc2,rho_bd2,rho_cd2;
    real drho_ab2,drho_ac2,drho_ad2,drho_bc2,drho_bd2,drho_cd2;
    real dval1,dval2,dval3,dval4,dval5,dval6;
    real dval1b,dval2b,dval3b,dval4b,dval5b,dval6b;
    real cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd;
    real cosine_abc[3],cosine_abd[3],cosine_acd[3],cosine_bcd[3];
    real invsin[6],cotan[6];
    real deriv_abc[3][3],deriv_abd[3][3];
    real deriv_acd[3][3],deriv_bcd[3][3];
    real dinvsin[6][6],dcotan[6][6];
    real dval1_ab[6],dval1_ac[6],dval1_ad[6],dval1_bc[6];
    real dval1_bd[6],dval1_cd[6];
    real dval2_ab[6],dval2_ac[6],dval2_ad[6],dval2_bc[6];
    real dval2_bd[6],dval2_cd[6];
    real dcap_ab[6],dcap_ac[6],dcap_ad[6],dcap_bc[6];
    real dcap_bd[6],dcap_cd[6];

    l1 = plane_dist(ra2, rb2, rab2);
    l2 = plane_dist(ra2, rc2, rac2);
    l3 = plane_dist(ra2, rd2, rad2);
    l4 = plane_dist(rb2, rc2, rbc2);
    l5 = plane_dist(rb2, rd2, rbd2);
    l6 = plane_dist(rc2, rd2, rcd2);

    val1 = l1*rab; val2 = l2*rac; val3 = l3*rad;
    val4 = l4*rbc; val5 = l5*rbd; val6 = l6*rcd;

    val1b = rab-val1; val2b = rac-val2; val3b = rad-val3;
    val4b = rbc-val4; val5b = rbd-val5; val6b = rcd-val6;

    // We consider the tetrahedron (A,B,C,P_ABC) where P_ABC is the
    // point of intersection of the three spheres such that (A,B,C,P_ABC) is ccw.
    // The edge lengths in this tetrahedron are: rab, rac, rAP=ra, rbc, rBP=rb, rCP=rc

    tet3dihedcos<compder>(rab2, rac2, ra2, rbc2, rb2, rc2, cosine_abc, deriv_abc);

    // repeat for tetrahedron (A,B,D,P_ABD)
    tet3dihedcos<compder>(rab2, rad2, ra2, rbd2, rb2, rd2, cosine_abd, deriv_abd);

    // repeat for tetrahedron (A,C,D,P_ACD)
    tet3dihedcos<compder>(rac2, rad2, ra2, rcd2, rc2, rd2, cosine_acd, deriv_acd);

    // repeat for tetrahedron (B,C,D,P_BCD)
    tet3dihedcos<compder>(rbc2, rbd2, rb2, rcd2, rc2, rd2, cosine_bcd, deriv_bcd);

    cos_abc = cosine_abc[0]; cos_acb = cosine_abc[1]; cos_bca = cosine_abc[2];
    cos_abd = cosine_abd[0]; cos_adb = cosine_abd[1]; cos_bda = cosine_abd[2];
    cos_acd = cosine_acd[0]; cos_adc = cosine_acd[1]; cos_cda = cosine_acd[2];
    cos_bcd = cosine_bcd[0]; cos_bdc = cosine_bcd[1]; cos_cdb = cosine_bcd[2];

    rho_ab2 = ra2 - val1b*val1b; rho_ac2 = ra2 - val2b*val2b;
    rho_ad2 = ra2 - val3b*val3b; rho_bc2 = rb2 - val4b*val4b;
    rho_bd2 = rb2 - val5b*val5b; rho_cd2 = rc2 - val6b*val6b;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(sin_ang[i]) < eps) {
            invsin[i] = 0;
            cotan[i] = 0;
        }
        else {
            invsin[i] = 1.0/sin_ang[i];
            cotan[i] = cos_ang[i]*invsin[i];
        }
    }

    val_ab = -(cos_abc*cos_abc+cos_abd*cos_abd)*cotan[0]
        + 2*cos_abc*cos_abd*invsin[0];
    val_ac = -(cos_acb*cos_acb+cos_acd*cos_acd)*cotan[1]
        + 2*cos_acb*cos_acd*invsin[1];
    val_ad = -(cos_adb*cos_adb+cos_adc*cos_adc)*cotan[2]
        + 2*cos_adb*cos_adc*invsin[2];
    val_bc = -(cos_bca*cos_bca+cos_bcd*cos_bcd)*cotan[3]
        + 2*cos_bca*cos_bcd*invsin[3];
    val_bd = -(cos_bda*cos_bda+cos_bdc*cos_bdc)*cotan[4]
        + 2*cos_bda*cos_bdc*invsin[4];
    val_cd = -(cos_cda*cos_cda+cos_cdb*cos_cdb)*cotan[5]
        + 2*cos_cda*cos_cdb*invsin[5];

    cap_ab = rho_ab2*val_ab; cap_ac = rho_ac2*val_ac;
    cap_ad = rho_ad2*val_ad; cap_bc = rho_bc2*val_bc;
    cap_bd = rho_bd2*val_bd; cap_cd = rho_cd2*val_cd;

    vola = (val1b*cap_ab+val2b*cap_ac+val3b*cap_ad)/6;
    volb = (val1*cap_ab+val4b*cap_bc+val5b*cap_bd)/6;
    volc = (val2*cap_ac+val4*cap_bc+val6b*cap_cd)/6;
    vold = (val3*cap_ad+val5*cap_bd+val6*cap_cd)/6;

    if constexpr (!compder) return;

    dval1b = l1; dval2b = l2; dval3b = l3;
    dval4b = l4; dval5b = l5; dval6b = l6;

    dval1 = 1-l1; dval2 = 1-l2; dval3 = 1-l3;
    dval4 = 1-l4; dval5 = 1-l5; dval6 = 1-l6;

    drho_ab2 = -2*dval1b*val1b; drho_ac2 = -2*dval2b*val2b;
    drho_ad2 = -2*dval3b*val3b; drho_bc2 = -2*dval4b*val4b;
    drho_bd2 = -2*dval5b*val5b; drho_cd2 = -2*dval6b*val6b;

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            dcotan[i][j] = -deriv[i][j]*(1.+cotan[i]*cotan[i]);
            dinvsin[i][j] = -deriv[i][j]*cotan[i]*invsin[i];
        }
    }

    val1_ab = cos_abc*cos_abc+cos_abd*cos_abd;
    val2_ab = 2*cos_abc*cos_abd;

    dval1_ab[0] = 2*(deriv_abc[0][0]*cos_abc+deriv_abd[0][0]*cos_abd);
    dval1_ab[1] = 2*deriv_abc[0][1]*cos_abc;
    dval1_ab[2] = 2*deriv_abd[0][1]*cos_abd;
    dval1_ab[3] = 2*deriv_abc[0][2]*cos_abc;
    dval1_ab[4] = 2*deriv_abd[0][2]*cos_abd;
    dval1_ab[5] = 0;

    dval2_ab[0] = 2*(deriv_abc[0][0]*cos_abd+deriv_abd[0][0]*cos_abc);
    dval2_ab[1] = 2*deriv_abc[0][1]*cos_abd;
    dval2_ab[2] = 2*deriv_abd[0][1]*cos_abc;
    dval2_ab[3] = 2*deriv_abc[0][2]*cos_abd;
    dval2_ab[4] = 2*deriv_abd[0][2]*cos_abc;
    dval2_ab[5] = 0;

    for (int i = 0; i < 6; i++) {
        dcap_ab[i] = - dval1_ab[i]*cotan[0] - val1_ab*dcotan[0][i]
            + dval2_ab[i]*invsin[0] + val2_ab*dinvsin[0][i];
        dcap_ab[i] = rho_ab2*dcap_ab[i];
    }
    dcap_ab[0] += drho_ab2*val_ab;

    val1_ac = cos_acb*cos_acb + cos_acd*cos_acd;
    val2_ac = 2*cos_acb*cos_acd;

    dval1_ac[0] = 2*deriv_abc[1][0]*cos_acb;
    dval1_ac[1] = 2*(deriv_abc[1][1]*cos_acb+deriv_acd[0][0]*cos_acd);
    dval1_ac[2] = 2*deriv_acd[0][1]*cos_acd;
    dval1_ac[3] = 2*deriv_abc[1][2]*cos_acb;
    dval1_ac[4] = 0;
    dval1_ac[5] = 2*deriv_acd[0][2]*cos_acd;

    dval2_ac[0] = 2*deriv_abc[1][0]*cos_acd;
    dval2_ac[1] = 2*(deriv_abc[1][1]*cos_acd+deriv_acd[0][0]*cos_acb);
    dval2_ac[2] = 2*deriv_acd[0][1]*cos_acb;
    dval2_ac[3] = 2*deriv_abc[1][2]*cos_acd;
    dval2_ac[4] = 0;
    dval2_ac[5] = 2*deriv_acd[0][2]*cos_acb;

    for (int i = 0; i < 6; i++) {
        dcap_ac[i] = -dval1_ac[i]*cotan[1] - val1_ac*dcotan[1][i]
            + dval2_ac[i]*invsin[1] + val2_ac*dinvsin[1][i];
        dcap_ac[i] = rho_ac2*dcap_ac[i];
    }
    dcap_ac[1] += drho_ac2*val_ac;

    val1_ad = cos_adb*cos_adb + cos_adc*cos_adc;
    val2_ad = 2*cos_adb*cos_adc;

    dval1_ad[0] = 2*deriv_abd[1][0]*cos_adb;
    dval1_ad[1] = 2*deriv_acd[1][0]*cos_adc;
    dval1_ad[2] = 2*(deriv_abd[1][1]*cos_adb+deriv_acd[1][1]*cos_adc);
    dval1_ad[3] = 0;
    dval1_ad[4] = 2*deriv_abd[1][2]*cos_adb;
    dval1_ad[5] = 2*deriv_acd[1][2]*cos_adc;

    dval2_ad[0] = 2*deriv_abd[1][0]*cos_adc;
    dval2_ad[1] = 2*deriv_acd[1][0]*cos_adb;
    dval2_ad[2] = 2*(deriv_abd[1][1]*cos_adc+deriv_acd[1][1]*cos_adb);
    dval2_ad[3] = 0;
    dval2_ad[4] = 2*deriv_abd[1][2]*cos_adc;
    dval2_ad[5] = 2*deriv_acd[1][2]*cos_adb;

    for (int i = 0; i < 6; i++) {
        dcap_ad[i] = -dval1_ad[i]*cotan[2] - val1_ad*dcotan[2][i]
            + dval2_ad[i]*invsin[2] + val2_ad*dinvsin[2][i];
        dcap_ad[i] = rho_ad2*dcap_ad[i];
    }
    dcap_ad[2] += drho_ad2*val_ad;

    val1_bc = cos_bca*cos_bca + cos_bcd*cos_bcd;
    val2_bc = 2*cos_bca*cos_bcd;

    dval1_bc[0] = 2*deriv_abc[2][0]*cos_bca;
    dval1_bc[1] = 2*deriv_abc[2][1]*cos_bca;
    dval1_bc[2] = 0;
    dval1_bc[3] = 2*(deriv_abc[2][2]*cos_bca+deriv_bcd[0][0]*cos_bcd);
    dval1_bc[4] = 2*deriv_bcd[0][1]*cos_bcd;
    dval1_bc[5] = 2*deriv_bcd[0][2]*cos_bcd;

    dval2_bc[0] = 2*deriv_abc[2][0]*cos_bcd;
    dval2_bc[1] = 2*deriv_abc[2][1]*cos_bcd;
    dval2_bc[2] = 0;
    dval2_bc[3] = 2*(deriv_abc[2][2]*cos_bcd+deriv_bcd[0][0]*cos_bca);
    dval2_bc[4] = 2*deriv_bcd[0][1]*cos_bca;
    dval2_bc[5] = 2*deriv_bcd[0][2]*cos_bca;

    for (int i = 0; i < 6; i++) {
        dcap_bc[i] = -dval1_bc[i]*cotan[3] - val1_bc*dcotan[3][i]
            + dval2_bc[i]*invsin[3] + val2_bc*dinvsin[3][i];
        dcap_bc[i] = rho_bc2*dcap_bc[i];
    }
    dcap_bc[3] += drho_bc2*val_bc;

    val1_bd = cos_bda*cos_bda + cos_bdc*cos_bdc;
    val2_bd = 2*cos_bda*cos_bdc;

    dval1_bd[0] = 2*deriv_abd[2][0]*cos_bda;
    dval1_bd[1] = 0;
    dval1_bd[2] = 2*deriv_abd[2][1]*cos_bda;
    dval1_bd[3] = 2*deriv_bcd[1][0]*cos_bdc;
    dval1_bd[4] = 2*(deriv_abd[2][2]*cos_bda+deriv_bcd[1][1]*cos_bdc);
    dval1_bd[5] = 2*deriv_bcd[1][2]*cos_bdc;

    dval2_bd[0] = 2*deriv_abd[2][0]*cos_bdc;
    dval2_bd[1] = 0;
    dval2_bd[2] = 2*deriv_abd[2][1]*cos_bdc;
    dval2_bd[3] = 2*deriv_bcd[1][0]*cos_bda;
    dval2_bd[4] = 2*(deriv_abd[2][2]*cos_bdc+deriv_bcd[1][1]*cos_bda);
    dval2_bd[5] = 2*deriv_bcd[1][2]*cos_bda;

    for (int i = 0; i < 6; i++) {
        dcap_bd[i] = -dval1_bd[i]*cotan[4] - val1_bd*dcotan[4][i]
            + dval2_bd[i]*invsin[4] + val2_bd*dinvsin[4][i];
        dcap_bd[i] = rho_bd2*dcap_bd[i];
    }
    dcap_bd[4] += drho_bd2*val_bd;

    val1_cd = cos_cda*cos_cda + cos_cdb*cos_cdb;
    val2_cd = 2*cos_cda*cos_cdb;

    dval1_cd[0] = 0;
    dval1_cd[1] = 2*deriv_acd[2][0]*cos_cda;
    dval1_cd[2] = 2*deriv_acd[2][1]*cos_cda;
    dval1_cd[3] = 2*deriv_bcd[2][0]*cos_cdb;
    dval1_cd[4] = 2*deriv_bcd[2][1]*cos_cdb;
    dval1_cd[5] = 2*(deriv_acd[2][2]*cos_cda+deriv_bcd[2][2]*cos_cdb);

    dval2_cd[0] = 0;
    dval2_cd[1] = 2*deriv_acd[2][0]*cos_cdb;
    dval2_cd[2] = 2*deriv_acd[2][1]*cos_cdb;
    dval2_cd[3] = 2*deriv_bcd[2][0]*cos_cda;
    dval2_cd[4] = 2*deriv_bcd[2][1]*cos_cda;
    dval2_cd[5] = 2*(deriv_acd[2][2]*cos_cdb+deriv_bcd[2][2]*cos_cda);

    for (int i = 0; i < 6; i++) {
        dcap_cd[i] = -dval1_cd[i]*cotan[5]-val1_cd*dcotan[5][i]
            +dval2_cd[i]*invsin[5]+val2_cd*dinvsin[5][i];
        dcap_cd[i] = rho_cd2*dcap_cd[i];
    }
    dcap_cd[5] += drho_cd2*val_cd;

    for (int i = 0; i < 6; i++) {
        dvola[i] = (val1b*dcap_ab[i]+val2b*dcap_ac[i]
            + val3b*dcap_ad[i])/6;
        dvolb[i] = (val1*dcap_ab[i]+val4b*dcap_bc[i]
            + val5b*dcap_bd[i])/6;
        dvolc[i] = (val2*dcap_ac[i]+val4*dcap_bc[i]
            + val6b*dcap_cd[i])/6;
        dvold[i] = (val3*dcap_ad[i]+val5*dcap_bd[i]
            + val6*dcap_cd[i])/6;
    }

    dvola[0] += dval1b*cap_ab/6;
    dvola[1] += dval2b*cap_ac/6;
    dvola[2] += dval3b*cap_ad/6;
    dvolb[0] += dval1*cap_ab/6;
    dvolb[3] += dval4b*cap_bc/6;
    dvolb[4] += dval5b*cap_bd/6;
    dvolc[1] += dval2*cap_ac/6;
    dvolc[3] += dval4*cap_bc/6;
    dvolc[5] += dval6b*cap_cd/6;
    dvold[2] += dval3*cap_ad/6;
    dvold[4] += dval5*cap_bd/6;
    dvold[5] += dval6*cap_cd/6;
}

template <bool compder>
inline real AlphaVol::trig_dradius(real a, real b, real c, real *der_r)
{
    real u = 4*a*b*c - (a+b+c-1)*(a+b+c-1);
    real sqr_u = REAL_SQRT(u);

    real v = (a-1)*(a-1) + (b-1)*(b-1) + (c-1)*(c-1) - (a-b)*(a-b) - (a-c)*(a-c) - (b-c)*(b-c);
    real sqr_v = REAL_SQRT(v);

    real r = 0.5 + 0.5* sqr_u/sqr_v;

    if constexpr (!compder) return r;

    real du_da = 4*b*c - 2*(a+b+c-1);
    real du_db = 4*a*c - 2*(a+b+c-1);
    real du_dc = 4*a*b - 2*(a+b+c-1);

    real dv_da = 2*(b+c-a-1);
    real dv_db = 2*(a+c-b-1);
    real dv_dc = 2*(a+b-c-1);

    der_r[0] = 0.5*(r-0.5)*(du_da/u -dv_da/v);
    der_r[1] = 0.5*(r-0.5)*(du_db/u -dv_db/v);
    der_r[2] = 0.5*(r-0.5)*(du_dc/u -dv_dc/v);

    return r;
}

///////////////////////////////////////////////////////
//                                                   //
//  trig_darea  --  spherical triangle surface area  //
//                                                   //
///////////////////////////////////////////////////////

// "trig_darea" omputes the surface area S of a
// spherical triangle and its derivatives

template <bool compder>
inline real AlphaVol::trig_darea(real a, real b, real c, real *der_S)
{
    real tol = 1.e-14;

    real u = 4*a*b*c - (a+b+c-1)*(a+b+c-1);
    real v = 4*a*b*c;
    if (REAL_ABS(u) < tol) u = 0.0;
    real w = REAL_SQRT(REAL_ABS(u));
    real t = REAL_SQRT(v);
    real wt = w/t;

    if (REAL_ABS(wt - 1) < eps) wt = 1;
    else if (REAL_ABS(wt + 1) < eps) wt = -1;

    real S = 2*REAL_ASIN(wt);

    if constexpr (!compder) return S;

    if (w>0) {
        der_S[0] = (b+c-a-1)/(a*w);
        der_S[1] = (a+c-b-1)/(b*w);
        der_S[2] = (a+b-c-1)/(c*w);
    }
    else {
        der_S[0] = 0;
        der_S[1] = 0;
        der_S[2] = 0;
    }

    return S;
}

////////////////////////////////////////////
//                                        //
//  sign  --  sign of spherical triangle  //
//                                        //
////////////////////////////////////////////

// "sign" defines if the surface area of a
// spherical triangle is positive or negative

inline real AlphaVol::sign(real a, real b, real c)
{
    real s = -1;
    if ( a + b <= 1 + c) {
        s = 1;
    }
    else {
        s = -1.;
    }

    return s;
}

///////////////////////////////////////////////////////////
//                                                       //
//  threesphgss  --  contribution to gaussian curvature  //
//                                                       //
///////////////////////////////////////////////////////////

// "threesphgss" computes the relative contribution of the
// 3 vertices of a triangle in the alpha complex to the
// Gaussian curvature

template <bool compder>
inline void AlphaVol::threesphgss(real ra, real rb, real rc, 
    real ra2, real rb2, real rc2, real rab, real rac, real rbc,
    real rab2, real rac2, real rbc2, real& areaA, real& areaB, 
    real& areaC, real darea[3][3])
{
    real cos_ab, cos_ac, cos_bc;
    real da_ab, db_bc, dc_ac;
    real a, b, c;
    real r;
    real sign_a, sign_b, sign_c;
    real S, Sa, Sb, Sc;
    real der_r[3], der_S[3], der_Sa[3], der_Sb[3], der_Sc[3];
    real dSa[3], dSb[3], dSc[3];
    real dSa_dab, dSa_dac, dSa_dbc;
    real dSb_dab, dSb_dac, dSb_dbc;
    real dSc_dab, dSc_dac, dSc_dbc;

    cos_ab = (ra2+rb2-rab2)/(2.0*ra*rb);
    a = 0.5*(cos_ab+1);

    cos_bc = (rb2+rc2-rbc2)/(2.0*rb*rc);
    b = 0.5*(cos_bc+1);

    cos_ac = (ra2+rc2-rac2)/(2.0*ra*rc);
    c = 0.5*(cos_ac+1);

    sign_a = sign(a, c, b);
    sign_b = sign(a, b, c);
    sign_c = sign(c, b, a);

    der_r[0] = 0; der_r[1] = 0; der_r[2] = 0;
    r = trig_dradius<compder>(a, b, c, der_r);

    der_S[0] = 0; der_S[1] = 0; der_S[2] = 0;
    der_Sa[0] = 0; der_Sa[1] = 0; der_Sa[2] = 0;
    der_Sb[0] = 0; der_Sb[1] = 0; der_Sb[2] = 0;
    der_Sc[0] = 0; der_Sc[1] = 0; der_Sc[2] = 0;

    S  = trig_darea<compder>(a, b, c, der_S);
    Sa = trig_darea<compder>(a, r, r, der_Sa);
    Sb = trig_darea<compder>(r, b, r, der_Sb);
    Sc = trig_darea<compder>(r, r, c, der_Sc);

    dSa[0] = der_Sa[0] + (der_Sa[1]+der_Sa[2])*der_r[0];
    dSa[1] = (der_Sa[1]+der_Sa[2])*der_r[1];
    dSa[2] = (der_Sa[1]+der_Sa[2])*der_r[2];

    dSb[0] = (der_Sb[0]+der_Sb[2])*der_r[0];
    dSb[1] = der_Sb[1] + (der_Sb[0]+der_Sb[2])*der_r[1];
    dSb[2] = (der_Sb[0]+der_Sb[2])*der_r[2];

    dSc[0] = (der_Sc[0]+der_Sc[1])*der_r[0];
    dSc[1] = (der_Sc[0]+der_Sc[1])*der_r[1];
    dSc[2] = der_Sc[2] + (der_Sc[0]+der_Sc[1])*der_r[2];

    if (Sa==0) {
        Sa = S - sign_a*Sb - sign_b*Sc;
        dSa[0] = der_S[0] - sign_a*dSb[0] - sign_b*dSc[0];
        dSa[1] = der_S[1] - sign_a*dSb[1] - sign_b*dSc[1];
        dSa[2] = der_S[2] - sign_a*dSb[2] - sign_b*dSc[2];
    }
    if (Sb==0) {
        Sb = S - sign_c*Sa - sign_b*Sc;
        dSb[0] = der_S[0] - sign_c*dSa[0] - sign_b*dSc[0];
        dSb[1] = der_S[1] - sign_c*dSa[1] - sign_b*dSc[1];
        dSb[2] = der_S[2] - sign_c*dSa[2] - sign_b*dSc[2];
    }
    if (Sc==0) {
        Sc = S - sign_c*Sa - sign_a*Sb;
        dSc[0] = der_S[0] - sign_c*dSa[0] - sign_a*dSb[0];
        dSc[1] = der_S[1] - sign_c*dSa[1] - sign_a*dSb[1];
        dSc[2] = der_S[2] - sign_c*dSa[2] - sign_a*dSb[2];
    }

    areaA = 0.5*(sign_c*Sa + sign_b*Sc);
    areaB = 0.5*(sign_a*Sb + sign_c*Sa);
    areaC = 0.5*(sign_b*Sc + sign_a*Sb);

    if constexpr (!compder) return;

    da_ab = -0.5*rab/(ra*rb);
    db_bc = -0.5*rbc/(rb*rc);
    dc_ac = -0.5*rac/(ra*rc);

    dSa_dab = dSa[0]*da_ab;
    dSa_dbc = dSa[1]*db_bc;
    dSa_dac = dSa[2]*dc_ac;

    dSb_dab = dSb[0]*da_ab;
    dSb_dbc = dSb[1]*db_bc;
    dSb_dac = dSb[2]*dc_ac;

    dSc_dab = dSc[0]*da_ab;
    dSc_dbc = dSc[1]*db_bc;
    dSc_dac = dSc[2]*dc_ac;

    darea[0][0] = 0.5*(sign_c*dSa_dab + sign_b*dSc_dab);
    darea[0][1] = 0.5*(sign_c*dSa_dac + sign_b*dSc_dac);
    darea[0][2] = 0.5*(sign_c*dSa_dbc + sign_b*dSc_dbc);

    darea[1][0] = 0.5*(sign_a*dSb_dab + sign_c*dSa_dab);
    darea[1][1] = 0.5*(sign_a*dSb_dac + sign_c*dSa_dac);
    darea[1][2] = 0.5*(sign_a*dSb_dbc + sign_c*dSa_dbc);

    darea[2][0] = 0.5*(sign_b*dSc_dab + sign_a*dSb_dab);
    darea[2][1] = 0.5*(sign_b*dSc_dac + sign_a*dSb_dac);
    darea[2][2] = 0.5*(sign_b*dSc_dbc + sign_a*dSb_dbc);
}

// explicit instatiation
template void AlphaVol::alphavol<true>(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::vector<Edge>& edges, std::vector<Face>& faces,
    real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
    real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord);
template void AlphaVol::alphavol<false>(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::vector<Edge>& edges, std::vector<Face>& faces,
    real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
    real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord);
}
