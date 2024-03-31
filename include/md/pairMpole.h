// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "damping.h"
#include "libfunc.h"
#include "mplpot.h"
#include <cmath>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  pairMpole  --  pairwise electrostatic calculation  //
//                                                     //
/////////////////////////////////////////////////////////

// "pairMpole" calculates the pairwise electrostatic energy 
// and/or gradient due to atomic multipole interactions

template <bool do_e, bool do_g, bool do_v>
inline void pairMpole(
    real r2, real xr, real yr, real zr, real mscale,
    real ci, real dix, real diy, real diz, real qixx, real qixy, real qixz,
    real qiyy, real qiyz, real qizz,
    real ck, real dkx, real dky, real dkz, real qkxx, real qkxy, real qkxz,
    real qkyy, real qkyz, real qkzz,
    real f,
    real& frcx, real& frcy, real& frcz,
    real& ttmxi, real& ttmyi, real& ttmzi,
    real& ttmxk, real& ttmyk, real& ttmzk,
    real& e,
    real& vxx, real& vxy, real& vxz,
    real& vyy, real& vyz, real& vzz)
{
    real r = REAL_SQRT(r2);
    real invr1 = 1/r;
    real rr2 = invr1 * invr1;

    real rr1,rr3,rr5,rr7,rr9,rr11;
    rr1 = mscale * f * invr1;
    rr3 = rr1 * rr2;
    rr5 = 3 * rr3 * rr2;
    rr7 = 5 * rr5 * rr2;
    rr9 = 7 * rr7 * rr2;
    if constexpr (do_g) rr11 = 9 * rr9 * rr2;

    real dir = dix * xr + diy * yr + diz * zr;
    real qix = qixx * xr + qixy * yr + qixz * zr;
    real qiy = qixy * xr + qiyy * yr + qiyz * zr;
    real qiz = qixz * xr + qiyz * yr + qizz * zr;
    real qir = qix * xr + qiy * yr + qiz * zr;
    real dkr = dkx * xr + dky * yr + dkz * zr;
    real qkx = qkxx * xr + qkxy * yr + qkxz * zr;
    real qky = qkxy * xr + qkyy * yr + qkyz * zr;
    real qkz = qkxz * xr + qkyz * yr + qkzz * zr;
    real qkr = qkx * xr + qky * yr + qkz * zr;
    real dik = dix * dkx + diy * dky + diz * dkz;
    real qik = qix * qkx + qiy * qky + qiz * qkz;
    real diqk = dix * qkx + diy * qky + diz * qkz;
    real dkqi = dkx * qix + dky * qiy + dkz * qiz;
    real qiqk = 2 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx
        + qiyy * qkyy + qizz * qkzz;

    real term1 = ci * ck;
    real term2 = ck * dir - ci * dkr + dik;
    real term3 = ci * qkr + ck * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
    real term4 = dir * qkr - dkr * qir - 4 * qik;
    real term5 = qir * qkr;

    if constexpr (do_e) {
        e = term1 * rr1 + term2 * rr3 + term3 * rr5 + term4 * rr7 + term5 * rr9;
    }

    if constexpr (do_g) {
        // gradient
        real qixk = qixx * qkx + qixy * qky + qixz * qkz;
        real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
        real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
        real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz;
        real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
        real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

        real diqkx = dix * qkxx + diy * qkxy + diz * qkxz;
        real diqky = dix * qkxy + diy * qkyy + diz * qkyz;
        real diqkz = dix * qkxz + diy * qkyz + diz * qkzz;
        real dkqix = dkx * qixx + dky * qixy + dkz * qixz;
        real dkqiy = dkx * qixy + dky * qiyy + dkz * qiyz;
        real dkqiz = dkx * qixz + dky * qiyz + dkz * qizz;

        real de = term1 * rr3 + term2 * rr5 + term3 * rr7 + term4 * rr9
            + term5 * rr11;

        term1 = -ck * rr3 + dkr * rr5 - qkr * rr7;
        term2 = ci * rr3 + dir * rr5 + qir * rr7;
        term3 = 2 * rr5;
        term4 = 2 * (-ck * rr5 + dkr * rr7 - qkr * rr9);
        term5 = 2 * (-ci * rr5 - dir * rr7 - qir * rr9);
        real term6 = 4 * rr7;

        frcx = de * xr + term1 * dix + term2 * dkx + term3 * (diqkx - dkqix)
            + term4 * qix + term5 * qkx + term6 * (qixk + qkxi);
        frcy = de * yr + term1 * diy + term2 * dky + term3 * (diqky - dkqiy)
            + term4 * qiy + term5 * qky + term6 * (qiyk + qkyi);
        frcz = de * zr + term1 * diz + term2 * dkz + term3 * (diqkz - dkqiz)
            + term4 * qiz + term5 * qkz + term6 * (qizk + qkzi);

        // torque
        real dirx = diy * zr - diz * yr;
        real diry = diz * xr - dix * zr;
        real dirz = dix * yr - diy * xr;
        real dkrx = dky * zr - dkz * yr;
        real dkry = dkz * xr - dkx * zr;
        real dkrz = dkx * yr - dky * xr;
        real dikx = diy * dkz - diz * dky;
        real diky = diz * dkx - dix * dkz;
        real dikz = dix * dky - diy * dkx;

        real qirx = qiz * yr - qiy * zr;
        real qiry = qix * zr - qiz * xr;
        real qirz = qiy * xr - qix * yr;
        real qkrx = qkz * yr - qky * zr;
        real qkry = qkx * zr - qkz * xr;
        real qkrz = qky * xr - qkx * yr;
        real qikx = qky * qiz - qkz * qiy;
        real qiky = qkz * qix - qkx * qiz;
        real qikz = qkx * qiy - qky * qix;

        real qikrx = qizk * yr - qiyk * zr;
        real qikry = qixk * zr - qizk * xr;
        real qikrz = qiyk * xr - qixk * yr;
        real qkirx = qkzi * yr - qkyi * zr;
        real qkiry = qkxi * zr - qkzi * xr;
        real qkirz = qkyi * xr - qkxi * yr;

        real diqkrx = diqkz * yr - diqky * zr;
        real diqkry = diqkx * zr - diqkz * xr;
        real diqkrz = diqky * xr - diqkx * yr;
        real dkqirx = dkqiz * yr - dkqiy * zr;
        real dkqiry = dkqix * zr - dkqiz * xr;
        real dkqirz = dkqiy * xr - dkqix * yr;

        real dqikx = diy * qkz - diz * qky + dky * qiz - dkz * qiy
            - 2 * (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy
            - qiyz * qkyy - qizz * qkyz);
        real dqiky = diz * qkx - dix * qkz + dkz * qix - dkx * qiz
            - 2 * (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz
            - qixy * qkyz - qixz * qkzz);
        real dqikz = dix * qky - diy * qkx + dkx * qiy - dky * qix
            - 2 * (qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx
            - qiyy * qkxy - qiyz * qkxz);

        ttmxi = -rr3 * dikx + term1 * dirx + term3 * (dqikx + dkqirx)
            - term4 * qirx - term6 * (qikrx + qikx);
        ttmyi = -rr3 * diky + term1 * diry + term3 * (dqiky + dkqiry)
            - term4 * qiry - term6 * (qikry + qiky);
        ttmzi = -rr3 * dikz + term1 * dirz + term3 * (dqikz + dkqirz)
            - term4 * qirz - term6 * (qikrz + qikz);
        ttmxk = rr3 * dikx + term2 * dkrx - term3 * (dqikx + diqkrx)
            - term5 * qkrx - term6 * (qkirx - qikx);
        ttmyk = rr3 * diky + term2 * dkry - term3 * (dqiky + diqkry)
            - term5 * qkry - term6 * (qkiry - qiky);
        ttmzk = rr3 * dikz + term2 * dkrz - term3 * (dqikz + diqkrz)
            - term5 * qkrz - term6 * (qkirz - qikz);

        if constexpr (do_v) {
            vxx = -xr * frcx;
            vxy = -0.5 * (yr * frcx + xr * frcy);
            vxz = -0.5 * (zr * frcx + xr * frcz);
            vyy = -yr * frcy;
            vyz = -0.5 * (zr * frcy + yr * frcz);
            vzz = -zr * frcz;
        }
    }
}

////////////////////////////////////////////////////////////////////
//                                                                //
//  pairMpoleC  --  pairwise electrostatic chargepen calculation  //
//                                                                //
////////////////////////////////////////////////////////////////////

// "pairMpoleCP" calculates the pairwise electrostatic charge penetration
// energy and/or gradient due to atomic multipole interactions

template <bool do_e, bool do_g, bool do_v, PenTyp PenType, bool use_cf>
inline void pairMpoleCP(
    real r2, real xr, real yr, real zr, real mscale,
    real corei, real vali, real alphai,
    real dix, real diy, real diz, real qixx, real qixy, real qixz,
    real qiyy, real qiyz, real qizz,
    real corek, real valk, real alphak,
    real dkx, real dky, real dkz, real qkxx, real qkxy, real qkxz,
    real qkyy, real qkyz, real qkzz,
    real f,
    real& frcx, real& frcy, real& frcz,
    real& ttmxi, real& ttmyi, real& ttmzi,
    real& ttmxk, real& ttmyk, real& ttmzk,
    real& e,
    real& vxx, real& vxy, real& vxz,
    real& vyy, real& vyz, real& vzz,
    real& poti, real& potk)
{
    real r = REAL_SQRT(r2);
    real invr1 = 1/r;
    real rr2 = invr1 * invr1;

    real rr1,rr3,rr5,rr7,rr9,rr11;
    rr1 = mscale * f * invr1;
    rr3 = rr1 * rr2;
    rr5 = 3 * rr3 * rr2;
    rr7 = 5 * rr5 * rr2;
    rr9 = 7 * rr7 * rr2;
    if constexpr (do_g) rr11 = 9 * rr9 * rr2;

    real dir = dix * xr + diy * yr + diz * zr;
    real qix = qixx * xr + qixy * yr + qixz * zr;
    real qiy = qixy * xr + qiyy * yr + qiyz * zr;
    real qiz = qixz * xr + qiyz * yr + qizz * zr;
    real qir = qix * xr + qiy * yr + qiz * zr;
    real dkr = dkx * xr + dky * yr + dkz * zr;
    real qkx = qkxx * xr + qkxy * yr + qkxz * zr;
    real qky = qkxy * xr + qkyy * yr + qkyz * zr;
    real qkz = qkxz * xr + qkyz * yr + qkzz * zr;
    real qkr = qkx * xr + qky * yr + qkz * zr;
    real dik = dix * dkx + diy * dky + diz * dkz;
    real qik = qix * qkx + qiy * qky + qiz * qkz;
    real diqk = dix * qkx + diy * qky + diz * qkz;
    real dkqi = dkx * qix + dky * qiy + dkz * qiz;
    real qiqk = 2 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx
        + qiyy * qkyy + qizz * qkzz;

    real term1 = corei * corek;
    real term1i = corek * vali;
    real term2i = corek * dir;
    real term3i = corek * qir;
    real term1k = corei * valk;
    real term2k = -corei * dkr;
    real term3k = corei * qkr;
    real term1ik = vali * valk;
    real term2ik = valk * dir - vali * dkr + dik;
    real term3ik = vali * qkr + valk * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
    real term4ik = dir * qkr - dkr * qir - 4 * qik;
    real term5ik = qir * qkr;

    real dmpi[5],dmpk[5],dmpik[6];
    if constexpr (do_g) {
        damppole<11,PenType>(r,alphai,alphak,dmpi,dmpk,dmpik);
    }
    else {
        damppole<9,PenType>(r,alphai,alphak,dmpi,dmpk,dmpik);
    }
    real rr1i = dmpi[0] * rr1;
    real rr3i = dmpi[1] * rr3;
    real rr5i = dmpi[2] * rr5;
    real rr1k = dmpk[0] * rr1;
    real rr3k = dmpk[1] * rr3;
    real rr5k = dmpk[2] * rr5;
    real rr1ik = dmpik[0] * rr1;
    real rr3ik = dmpik[1] * rr3;
    real rr5ik = dmpik[2] * rr5;
    real rr7ik = dmpik[3] * rr7;
    real rr9ik = dmpik[4] * rr9;
    real rr7i,rr7k,rr11ik;
    if constexpr (do_g) {
        rr7i = dmpi[3] * rr7;
        rr7k = dmpk[3] * rr7;
        rr11ik = dmpik[5] * rr11;
    }

    if constexpr (do_e) {
        e = term1 * rr1 + term4ik * rr7ik + term5ik * rr9ik + term1i * rr1i
            + term1k * rr1k + term1ik * rr1ik + term2i * rr3i + term2k * rr3k
            + term2ik * rr3ik + term3i * rr5i + term3k * rr5k + term3ik * rr5ik;
    }

    if constexpr (do_g) {
        // gradient
        real qixk = qixx * qkx + qixy * qky + qixz * qkz;
        real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
        real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
        real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz;
        real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
        real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

        real diqkx = dix * qkxx + diy * qkxy + diz * qkxz;
        real diqky = dix * qkxy + diy * qkyy + diz * qkyz;
        real diqkz = dix * qkxz + diy * qkyz + diz * qkzz;
        real dkqix = dkx * qixx + dky * qixy + dkz * qixz;
        real dkqiy = dkx * qixy + dky * qiyy + dkz * qiyz;
        real dkqiz = dkx * qixz + dky * qiyz + dkz * qizz;

        real de = term1 * rr3 + term4ik * rr9ik + term5ik * rr11ik + term1i * rr3i
            + term1k * rr3k + term1ik * rr3ik + term2i * rr5i + term2k * rr5k
            + term2ik * rr5ik + term3i * rr7i + term3k * rr7k + term3ik * rr7ik;

        term1 = -corek * rr3i - valk * rr3ik + dkr * rr5ik - qkr * rr7ik;
        real term2 = corei * rr3k + vali * rr3ik + dir * rr5ik + qir * rr7ik;
        real term3 = 2 * rr5ik;
        real term4 = -2 * (corek * rr5i + valk * rr5ik - dkr * rr7ik + qkr * rr9ik);
        real term5 = -2 * (corei * rr5k + vali * rr5ik + dir * rr7ik + qir * rr9ik);
        real term6 = 4 * rr7ik;

        if constexpr (use_cf) {
            real t1i = corek * rr1i + valk * rr1ik;
            real t1k = corei * rr1k + vali * rr1ik;
            real t2i = -dkr * rr3ik;
            real t2k = dir * rr3ik;
            real t3i = qkr * rr5ik;
            real t3k = qir * rr5ik;
            poti = t1i + t2i + t3i;
            potk = t1k + t2k + t3k;
        }

        frcx = de * xr + term1 * dix + term2 * dkx + term3 * (diqkx - dkqix)
            + term4 * qix + term5 * qkx + term6 * (qixk + qkxi);
        frcy = de * yr + term1 * diy + term2 * dky + term3 * (diqky - dkqiy)
            + term4 * qiy + term5 * qky + term6 * (qiyk + qkyi);
        frcz = de * zr + term1 * diz + term2 * dkz + term3 * (diqkz - dkqiz)
            + term4 * qiz + term5 * qkz + term6 * (qizk + qkzi);

        // torque
        real dirx = diy * zr - diz * yr;
        real diry = diz * xr - dix * zr;
        real dirz = dix * yr - diy * xr;
        real dkrx = dky * zr - dkz * yr;
        real dkry = dkz * xr - dkx * zr;
        real dkrz = dkx * yr - dky * xr;
        real dikx = diy * dkz - diz * dky;
        real diky = diz * dkx - dix * dkz;
        real dikz = dix * dky - diy * dkx;

        real qirx = qiz * yr - qiy * zr;
        real qiry = qix * zr - qiz * xr;
        real qirz = qiy * xr - qix * yr;
        real qkrx = qkz * yr - qky * zr;
        real qkry = qkx * zr - qkz * xr;
        real qkrz = qky * xr - qkx * yr;
        real qikx = qky * qiz - qkz * qiy;
        real qiky = qkz * qix - qkx * qiz;
        real qikz = qkx * qiy - qky * qix;

        real qikrx = qizk * yr - qiyk * zr;
        real qikry = qixk * zr - qizk * xr;
        real qikrz = qiyk * xr - qixk * yr;
        real qkirx = qkzi * yr - qkyi * zr;
        real qkiry = qkxi * zr - qkzi * xr;
        real qkirz = qkyi * xr - qkxi * yr;

        real diqkrx = diqkz * yr - diqky * zr;
        real diqkry = diqkx * zr - diqkz * xr;
        real diqkrz = diqky * xr - diqkx * yr;
        real dkqirx = dkqiz * yr - dkqiy * zr;
        real dkqiry = dkqix * zr - dkqiz * xr;
        real dkqirz = dkqiy * xr - dkqix * yr;

        real dqikx = diy * qkz - diz * qky + dky * qiz - dkz * qiy
            - 2 * (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy
            - qiyz * qkyy - qizz * qkyz);
        real dqiky = diz * qkx - dix * qkz + dkz * qix - dkx * qiz
            - 2 * (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz
            - qixy * qkyz - qixz * qkzz);
        real dqikz = dix * qky - diy * qkx + dkx * qiy - dky * qix
            - 2 * (qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx
            - qiyy * qkxy - qiyz * qkxz);

        ttmxi = -rr3ik * dikx + term1 * dirx + term3 * (dqikx + dkqirx)
            - term4 * qirx - term6 * (qikrx + qikx);
        ttmyi = -rr3ik * diky + term1 * diry + term3 * (dqiky + dkqiry)
            - term4 * qiry - term6 * (qikry + qiky);
        ttmzi = -rr3ik * dikz + term1 * dirz + term3 * (dqikz + dkqirz)
            - term4 * qirz - term6 * (qikrz + qikz);
        ttmxk = rr3ik * dikx + term2 * dkrx - term3 * (dqikx + diqkrx)
            - term5 * qkrx - term6 * (qkirx - qikx);
        ttmyk = rr3ik * diky + term2 * dkry - term3 * (dqiky + diqkry)
            - term5 * qkry - term6 * (qkiry - qiky);
        ttmzk = rr3ik * dikz + term2 * dkrz - term3 * (dqikz + diqkrz)
            - term5 * qkrz - term6 * (qkirz - qikz);

        if constexpr (do_v) {
            vxx = -xr * frcx;
            vxy = -0.5 * (yr * frcx + xr * frcy);
            vxz = -0.5 * (zr * frcx + xr * frcz);
            vyy = -yr * frcy;
            vyz = -0.5 * (zr * frcy + yr * frcz);
            vzz = -zr * frcz;
        }
    }
}
}
