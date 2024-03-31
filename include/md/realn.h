// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
struct float3
{
   float x, y, z;
};

inline float3 make_float3(float x, float y, float z)
{
   float3 f{x, y, z};
   return f;
}

struct double3
{
   double x, y, z;
};

inline double3 make_double3(double x, double y, double z)
{
   double3 f{x, y, z};
   return f;
}

#if POLMDQC_REAL_SIZE == 8
using real3 = double3;
#   define make_real3(x, y, z)    make_double3((x), (y), (z))
#elif POLMDQC_REAL_SIZE == 4
using real3 = float3;
#   define make_real3(x, y, z)    make_float3((x), (y), (z))
#endif

inline real3 operator-(real3 a)
{
   return make_real3(-a.x, -a.y, -a.z);
}

inline real3 operator+(real3 a, real3 b)
{
   return make_real3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline real3 operator-(real3 a, real3 b)
{
   return make_real3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline real3 operator*(real a, real3 b)
{
   return make_real3(a * b.x, a * b.y, a * b.z);
}

inline real3 operator*(real3 b, real a)
{
   return make_real3(a * b.x, a * b.y, a * b.z);
}

inline void operator+=(real3& a, real3 b)
{
   a.x += b.x;
   a.y += b.y;
   a.z += b.z;
}

inline void operator-=(real3& a, real3 b)
{
   a.x -= b.x;
   a.y -= b.y;
   a.z -= b.z;
}

inline void operator*=(real3& a, real b)
{
   a.x *= b;
   a.y *= b;
   a.z *= b;
}

inline real dot3(real3 a, real3 b)
{
   return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline real dot3(real3 a, real bx, real by, real bz)
{
   return a.x * bx + a.y * by + a.z * bz;
}

inline real dot3(real ax, real ay, real az, real3 b)
{
   return ax * b.x + ay * b.y + az * b.z;
}

inline real dot3(real ax, real ay, real az, real bx, real by, real bz)
{
   return ax * bx + ay * by + az * bz;
}

inline real3 cross(real ax, real ay, real az, real bx, real by, real bz)
{
   return make_real3(ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx);
}

inline real3 cross(real3 a, real3 b)
{
   return cross(a.x, a.y, a.z, b.x, b.y, b.z);
}

inline real3 cross(real3 a, real bx, real by, real bz)
{
   return cross(a.x, a.y, a.z, bx, by, bz);
}

inline real3 cross(real ax, real ay, real az, real3 b)
{
   return cross(ax, ay, az, b.x, b.y, b.z);
}
}
