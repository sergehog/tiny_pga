// 3D Projective Geometric Algebra
// Written by a generator written by enki.

/*
 * This file is part of the Tiny-PGA distribution (https://github.com/sergehog/tiny_pga)
 * Copyright (c) 2020 Sergey Smirnov / Seregium Oy.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <array>
#include <cmath>

#define PI 3.14159265358979323846

enum Basis : std::size_t
{
    kScalar = 0U,
    kE0 = 1U,
    kE1 = 2U,
    kNX = kE1,  // alias for plane x-normal
    kE2 = 3U,
    kNY = kE2,  // alias for plane y-normal
    kE3 = 4U,
    kNZ = kE3,  // alias for plane z-normal
    kE01 = 5U,  // X offset for line
    kE02 = 6U,  // Y offset for line
    kE03 = 7U,  // Z offset for line
    kE12 = 8U,
    kAX = kE12,  // bivector for X direction ?? WTF
    kE31 = 9U,
    kAY = kE12,  // bivector for Y direction
    kE23 = 10U,
    kAZ = kE23,  // bivector for Z direction ?? WTF
    kE021 = 11U,
    kX = kE23,  // X coordinate for a point
    kE013 = 12U,
    kY = kE23,  // Y coordinate for a point
    kE032 = 13U,
    kZ = kE23,  // Z coordinate for a point
    kE123 = 14U,
    kE0123 = 15U
};

template <typename ScalarType = float>
class PGA3D
{
  private:
    std::array<ScalarType, 16> mvec;

  public:
    PGA3D() : mvec{} {}

    PGA3D(ScalarType f, Basis idx = kScalar) : mvec{} { mvec[idx] = f; }

    PGA3D(Basis idx) : mvec{} { mvec[idx] = 1.F; }

    ScalarType& operator[](size_t idx) { return mvec[idx]; }
    const ScalarType& operator[](size_t idx) const { return mvec[idx]; }

    /// Clifford Conjugation  : res = a.Conjugate()
    PGA3D<ScalarType> Conjugate()
    {
        PGA3D<ScalarType> res;
        res[0] = mvec[0];
        res[1] = -mvec[1];
        res[2] = -mvec[2];
        res[3] = -mvec[3];
        res[4] = -mvec[4];
        res[5] = -mvec[5];
        res[6] = -mvec[6];
        res[7] = -mvec[7];
        res[8] = -mvec[8];
        res[9] = -mvec[9];
        res[10] = -mvec[10];
        res[11] = mvec[11];
        res[12] = mvec[12];
        res[13] = mvec[13];
        res[14] = mvec[14];
        res[15] = mvec[15];
        return res;
    };

    /// Main involution : res = a.Involute()
    PGA3D<ScalarType> Involute()
    {
        PGA3D<ScalarType> res;
        res[0] = mvec[0];
        res[1] = -mvec[1];
        res[2] = -mvec[2];
        res[3] = -mvec[3];
        res[4] = -mvec[4];
        res[5] = mvec[5];
        res[6] = mvec[6];
        res[7] = mvec[7];
        res[8] = mvec[8];
        res[9] = mvec[9];
        res[10] = mvec[10];
        res[11] = -mvec[11];
        res[12] = -mvec[12];
        res[13] = -mvec[13];
        res[14] = -mvec[14];
        res[15] = mvec[15];
        return res;
    };

    PGA3D<ScalarType> sqrt() { return (1.F + *this).normalized(); }

    float norm() { return std::sqrt(std::abs(float(((*this) * Conjugate()).mvec[0]))); }

    float inorm() { return (!(*this)).norm(); }

    PGA3D<ScalarType> normalized() { return (*this) * ScalarType(1.F / norm()); }
};

/// Reverse the order of the basis blades. : res = ~a
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator~(const PGA3D<ScalarType>& a)
{
    PGA3D<ScalarType> res;
    res[0] = a[0];
    res[1] = a[1];
    res[2] = a[2];
    res[3] = a[3];
    res[4] = a[4];
    res[5] = -a[5];
    res[6] = -a[6];
    res[7] = -a[7];
    res[8] = -a[8];
    res[9] = -a[9];
    res[10] = -a[10];
    res[11] = -a[11];
    res[12] = -a[12];
    res[13] = -a[13];
    res[14] = -a[14];
    res[15] = a[15];
    return res;
};

/// Poincare duality operator: res = !a
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator!(const PGA3D<ScalarType>& a)
{
    PGA3D<ScalarType> res;
    res[0] = a[15];
    res[1] = a[14];
    res[2] = a[13];
    res[3] = a[12];
    res[4] = a[11];
    res[5] = a[10];
    res[6] = a[9];
    res[7] = a[8];
    res[8] = a[7];
    res[9] = a[6];
    res[10] = a[5];
    res[11] = a[4];
    res[12] = a[3];
    res[13] = a[2];
    res[14] = a[1];
    res[15] = a[0];
    return res;
};

/// The geometric product: res = a * b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator*(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res{};
    res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] - b[9] * a[9] - b[10] * a[10] -
             b[14] * a[14];
    res[1] = b[1] * a[0] + b[0] * a[1] + b[2] * a[5] - b[5] * a[2] + b[3] * a[6] - b[6] * a[3] + b[4] * a[7] -
             b[7] * a[4] + b[8] * a[11] + b[11] * a[8] + b[9] * a[12] + b[12] * a[9] + b[10] * a[13] + b[13] * a[10] -
             b[14] * a[15] + b[15] * a[14];
    res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] - b[4] * a[9] - b[14] * a[10] -
             b[10] * a[14];
    res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] - b[2] * a[8] - b[14] * a[9] + b[4] * a[10] -
             b[9] * a[14];
    res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] - b[14] * a[8] + b[2] * a[9] - b[3] * a[10] -
             b[8] * a[14];
    res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] - b[11] * a[3] + b[12] * a[4] + b[0] * a[5] - b[8] * a[6] +
             b[9] * a[7] + b[6] * a[8] - b[7] * a[9] - b[15] * a[10] - b[3] * a[11] + b[4] * a[12] + b[14] * a[13] -
             b[13] * a[14] - b[10] * a[15];
    res[6] = b[6] * a[0] + b[3] * a[1] + b[11] * a[2] - b[1] * a[3] - b[13] * a[4] + b[8] * a[5] + b[0] * a[6] -
             b[10] * a[7] - b[5] * a[8] - b[15] * a[9] + b[7] * a[10] + b[2] * a[11] + b[14] * a[12] - b[4] * a[13] -
             b[12] * a[14] - b[9] * a[15];
    res[7] = b[7] * a[0] + b[4] * a[1] - b[12] * a[2] + b[13] * a[3] - b[1] * a[4] - b[9] * a[5] + b[10] * a[6] +
             b[0] * a[7] - b[15] * a[8] + b[5] * a[9] - b[6] * a[10] + b[14] * a[11] - b[2] * a[12] + b[3] * a[13] -
             b[11] * a[14] - b[8] * a[15];
    res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[14] * a[4] + b[0] * a[8] + b[10] * a[9] - b[9] * a[10] +
             b[4] * a[14];
    res[9] = b[9] * a[0] - b[4] * a[2] + b[14] * a[3] + b[2] * a[4] - b[10] * a[8] + b[0] * a[9] + b[8] * a[10] +
             b[3] * a[14];
    res[10] = b[10] * a[0] + b[14] * a[2] + b[4] * a[3] - b[3] * a[4] + b[9] * a[8] - b[8] * a[9] + b[0] * a[10] +
              b[2] * a[14];
    res[11] = b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] + b[15] * a[4] - b[3] * a[5] + b[2] * a[6] -
              b[14] * a[7] - b[1] * a[8] + b[13] * a[9] - b[12] * a[10] + b[0] * a[11] + b[10] * a[12] - b[9] * a[13] +
              b[7] * a[14] - b[4] * a[15];
    res[12] = b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[15] * a[3] + b[5] * a[4] + b[4] * a[5] - b[14] * a[6] -
              b[2] * a[7] - b[13] * a[8] - b[1] * a[9] + b[11] * a[10] - b[10] * a[11] + b[0] * a[12] + b[8] * a[13] +
              b[6] * a[14] - b[3] * a[15];
    res[13] = b[13] * a[0] - b[10] * a[1] + b[15] * a[2] + b[7] * a[3] - b[6] * a[4] - b[14] * a[5] - b[4] * a[6] +
              b[3] * a[7] + b[12] * a[8] - b[11] * a[9] - b[1] * a[10] + b[9] * a[11] - b[8] * a[12] + b[0] * a[13] +
              b[5] * a[14] - b[2] * a[15];
    res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] + b[4] * a[8] + b[3] * a[9] + b[2] * a[10] +
              b[0] * a[14];
    res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] + b[11] * a[4] + b[10] * a[5] + b[9] * a[6] +
              b[8] * a[7] + b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] - b[3] * a[12] - b[2] * a[13] -
              b[1] * a[14] + b[0] * a[15];
    return res;
};

/// The outer  product. Also known as Wedge or MEET : res = a ^ b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator^(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = b[0] * a[0];
    res[1] = b[1] * a[0] + b[0] * a[1];
    res[2] = b[2] * a[0] + b[0] * a[2];
    res[3] = b[3] * a[0] + b[0] * a[3];
    res[4] = b[4] * a[0] + b[0] * a[4];
    res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] + b[0] * a[5];
    res[6] = b[6] * a[0] + b[3] * a[1] - b[1] * a[3] + b[0] * a[6];
    res[7] = b[7] * a[0] + b[4] * a[1] - b[1] * a[4] + b[0] * a[7];
    res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[0] * a[8];
    res[9] = b[9] * a[0] - b[4] * a[2] + b[2] * a[4] + b[0] * a[9];
    res[10] = b[10] * a[0] + b[4] * a[3] - b[3] * a[4] + b[0] * a[10];
    res[11] =
        b[11] * a[0] - b[8] * a[1] + b[6] * a[2] - b[5] * a[3] - b[3] * a[5] + b[2] * a[6] - b[1] * a[8] + b[0] * a[11];
    res[12] =
        b[12] * a[0] - b[9] * a[1] - b[7] * a[2] + b[5] * a[4] + b[4] * a[5] - b[2] * a[7] - b[1] * a[9] + b[0] * a[12];
    res[13] = b[13] * a[0] - b[10] * a[1] + b[7] * a[3] - b[6] * a[4] - b[4] * a[6] + b[3] * a[7] - b[1] * a[10] +
              b[0] * a[13];
    res[14] = b[14] * a[0] + b[10] * a[2] + b[9] * a[3] + b[8] * a[4] + b[4] * a[8] + b[3] * a[9] + b[2] * a[10] +
              b[0] * a[14];
    res[15] = b[15] * a[0] + b[14] * a[1] + b[13] * a[2] + b[12] * a[3] + b[11] * a[4] + b[10] * a[5] + b[9] * a[6] +
              b[8] * a[7] + b[7] * a[8] + b[6] * a[9] + b[5] * a[10] - b[4] * a[11] - b[3] * a[12] - b[2] * a[13] -
              b[1] * a[14] + b[0] * a[15];
    return res;
};

/// The regressive product. Also known as Vee or JOIN : res = a & b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator&(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[15] = (a[15] * b[15]);
    res[14] = -(a[14] * -1 * b[15] + a[15] * b[14] * -1);
    res[13] = -(a[13] * -1 * b[15] + a[15] * b[13] * -1);
    res[12] = -(a[12] * -1 * b[15] + a[15] * b[12] * -1);
    res[11] = -(a[11] * -1 * b[15] + a[15] * b[11] * -1);
    res[10] = (a[10] * b[15] + a[13] * -1 * b[14] * -1 - a[14] * -1 * b[13] * -1 + a[15] * b[10]);
    res[9] = (a[9] * b[15] + a[12] * -1 * b[14] * -1 - a[14] * -1 * b[12] * -1 + a[15] * b[9]);
    res[8] = (a[8] * b[15] + a[11] * -1 * b[14] * -1 - a[14] * -1 * b[11] * -1 + a[15] * b[8]);
    res[7] = (a[7] * b[15] + a[12] * -1 * b[13] * -1 - a[13] * -1 * b[12] * -1 + a[15] * b[7]);
    res[6] = (a[6] * b[15] - a[11] * -1 * b[13] * -1 + a[13] * -1 * b[11] * -1 + a[15] * b[6]);
    res[5] = (a[5] * b[15] + a[11] * -1 * b[12] * -1 - a[12] * -1 * b[11] * -1 + a[15] * b[5]);
    res[4] = (a[4] * b[15] - a[7] * b[14] * -1 + a[9] * b[13] * -1 - a[10] * b[12] * -1 - a[12] * -1 * b[10] +
              a[13] * -1 * b[9] - a[14] * -1 * b[7] + a[15] * b[4]);
    res[3] = (a[3] * b[15] - a[6] * b[14] * -1 - a[8] * b[13] * -1 + a[10] * b[11] * -1 + a[11] * -1 * b[10] -
              a[13] * -1 * b[8] - a[14] * -1 * b[6] + a[15] * b[3]);
    res[2] = (a[2] * b[15] - a[5] * b[14] * -1 + a[8] * b[12] * -1 - a[9] * b[11] * -1 - a[11] * -1 * b[9] +
              a[12] * -1 * b[8] - a[14] * -1 * b[5] + a[15] * b[2]);
    res[1] = (a[1] * b[15] + a[5] * b[13] * -1 + a[6] * b[12] * -1 + a[7] * b[11] * -1 + a[11] * -1 * b[7] +
              a[12] * -1 * b[6] + a[13] * -1 * b[5] + a[15] * b[1]);
    res[0] = (a[0] * b[15] + a[1] * b[14] * -1 + a[2] * b[13] * -1 + a[3] * b[12] * -1 + a[4] * b[11] * -1 +
              a[5] * b[10] + a[6] * b[9] + a[7] * b[8] + a[8] * b[7] + a[9] * b[6] + a[10] * b[5] - a[11] * -1 * b[4] -
              a[12] * -1 * b[3] - a[13] * -1 * b[2] - a[14] * -1 * b[1] + a[15] * b[0]);
    return res;
};

/// The inner product. (Dot) : res = a | b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator|(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = b[0] * a[0] + b[2] * a[2] + b[3] * a[3] + b[4] * a[4] - b[8] * a[8] - b[9] * a[9] - b[10] * a[10] -
             b[14] * a[14];
    res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] - b[7] * a[4] + b[2] * a[5] + b[3] * a[6] +
             b[4] * a[7] + b[11] * a[8] + b[12] * a[9] + b[13] * a[10] + b[8] * a[11] + b[9] * a[12] + b[10] * a[13] +
             b[15] * a[14] - b[14] * a[15];
    res[2] = b[2] * a[0] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] + b[3] * a[8] - b[4] * a[9] - b[14] * a[10] -
             b[10] * a[14];
    res[3] = b[3] * a[0] + b[8] * a[2] + b[0] * a[3] - b[10] * a[4] - b[2] * a[8] - b[14] * a[9] + b[4] * a[10] -
             b[9] * a[14];
    res[4] = b[4] * a[0] - b[9] * a[2] + b[10] * a[3] + b[0] * a[4] - b[14] * a[8] + b[2] * a[9] - b[3] * a[10] -
             b[8] * a[14];
    res[5] = b[5] * a[0] - b[11] * a[3] + b[12] * a[4] + b[0] * a[5] - b[15] * a[10] - b[3] * a[11] + b[4] * a[12] -
             b[10] * a[15];
    res[6] = b[6] * a[0] + b[11] * a[2] - b[13] * a[4] + b[0] * a[6] - b[15] * a[9] + b[2] * a[11] - b[4] * a[13] -
             b[9] * a[15];
    res[7] = b[7] * a[0] - b[12] * a[2] + b[13] * a[3] + b[0] * a[7] - b[15] * a[8] - b[2] * a[12] + b[3] * a[13] -
             b[8] * a[15];
    res[8] = b[8] * a[0] + b[14] * a[4] + b[0] * a[8] + b[4] * a[14];
    res[9] = b[9] * a[0] + b[14] * a[3] + b[0] * a[9] + b[3] * a[14];
    res[10] = b[10] * a[0] + b[14] * a[2] + b[0] * a[10] + b[2] * a[14];
    res[11] = b[11] * a[0] + b[15] * a[4] + b[0] * a[11] - b[4] * a[15];
    res[12] = b[12] * a[0] + b[15] * a[3] + b[0] * a[12] - b[3] * a[15];
    res[13] = b[13] * a[0] + b[15] * a[2] + b[0] * a[13] - b[2] * a[15];
    res[14] = b[14] * a[0] + b[0] * a[14];
    res[15] = b[15] * a[0] + b[0] * a[15];
    return res;
};

/// Multivector addition : res = a + b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator+(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = a[0] + b[0];
    res[1] = a[1] + b[1];
    res[2] = a[2] + b[2];
    res[3] = a[3] + b[3];
    res[4] = a[4] + b[4];
    res[5] = a[5] + b[5];
    res[6] = a[6] + b[6];
    res[7] = a[7] + b[7];
    res[8] = a[8] + b[8];
    res[9] = a[9] + b[9];
    res[10] = a[10] + b[10];
    res[11] = a[11] + b[11];
    res[12] = a[12] + b[12];
    res[13] = a[13] + b[13];
    res[14] = a[14] + b[14];
    res[15] = a[15] + b[15];
    return res;
};

/// Multivector subtraction : res = a - b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator-(const PGA3D<ScalarType>& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = a[0] - b[0];
    res[1] = a[1] - b[1];
    res[2] = a[2] - b[2];
    res[3] = a[3] - b[3];
    res[4] = a[4] - b[4];
    res[5] = a[5] - b[5];
    res[6] = a[6] - b[6];
    res[7] = a[7] - b[7];
    res[8] = a[8] - b[8];
    res[9] = a[9] - b[9];
    res[10] = a[10] - b[10];
    res[11] = a[11] - b[11];
    res[12] = a[12] - b[12];
    res[13] = a[13] - b[13];
    res[14] = a[14] - b[14];
    res[15] = a[15] - b[15];
    return res;
};

/// scalar/multivector multiplication : res = a * b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator*(const ScalarType& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = a * b[0];
    res[1] = a * b[1];
    res[2] = a * b[2];
    res[3] = a * b[3];
    res[4] = a * b[4];
    res[5] = a * b[5];
    res[6] = a * b[6];
    res[7] = a * b[7];
    res[8] = a * b[8];
    res[9] = a * b[9];
    res[10] = a * b[10];
    res[11] = a * b[11];
    res[12] = a * b[12];
    res[13] = a * b[13];
    res[14] = a * b[14];
    res[15] = a * b[15];
    return res;
};

/// multivector/scalar multiplication : res = a * b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator*(const PGA3D<ScalarType>& a, const ScalarType& b)
{
    PGA3D<ScalarType> res;
    res[0] = a[0] * b;
    res[1] = a[1] * b;
    res[2] = a[2] * b;
    res[3] = a[3] * b;
    res[4] = a[4] * b;
    res[5] = a[5] * b;
    res[6] = a[6] * b;
    res[7] = a[7] * b;
    res[8] = a[8] * b;
    res[9] = a[9] * b;
    res[10] = a[10] * b;
    res[11] = a[11] * b;
    res[12] = a[12] * b;
    res[13] = a[13] * b;
    res[14] = a[14] * b;
    res[15] = a[15] * b;
    return res;
};

/// scalar/multivector addition : res = a + b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator+(const float& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = a + b[0];
    res[1] = b[1];
    res[2] = b[2];
    res[3] = b[3];
    res[4] = b[4];
    res[5] = b[5];
    res[6] = b[6];
    res[7] = b[7];
    res[8] = b[8];
    res[9] = b[9];
    res[10] = b[10];
    res[11] = b[11];
    res[12] = b[12];
    res[13] = b[13];
    res[14] = b[14];
    res[15] = b[15];
    return res;
};

/// multivector/scalar addition : res = a + b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator+(const PGA3D<ScalarType>& a, const float& b)
{
    PGA3D<ScalarType> res;
    res[0] = a[0] + b;
    res[1] = a[1];
    res[2] = a[2];
    res[3] = a[3];
    res[4] = a[4];
    res[5] = a[5];
    res[6] = a[6];
    res[7] = a[7];
    res[8] = a[8];
    res[9] = a[9];
    res[10] = a[10];
    res[11] = a[11];
    res[12] = a[12];
    res[13] = a[13];
    res[14] = a[14];
    res[15] = a[15];
    return res;
};

/// scalar/multivector subtraction : res = a - b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator-(const float& a, const PGA3D<ScalarType>& b)
{
    PGA3D<ScalarType> res;
    res[0] = a - b[0];
    res[1] = -b[1];
    res[2] = -b[2];
    res[3] = -b[3];
    res[4] = -b[4];
    res[5] = -b[5];
    res[6] = -b[6];
    res[7] = -b[7];
    res[8] = -b[8];
    res[9] = -b[9];
    res[10] = -b[10];
    res[11] = -b[11];
    res[12] = -b[12];
    res[13] = -b[13];
    res[14] = -b[14];
    res[15] = -b[15];
    return res;
};

/// multivector/scalar subtraction : res = a - b
template <typename ScalarType = float>
inline PGA3D<ScalarType> operator-(const PGA3D<ScalarType>& a, const float& b)
{
    PGA3D<ScalarType> res;
    res[0] = a[0] - b;
    res[1] = a[1];
    res[2] = a[2];
    res[3] = a[3];
    res[4] = a[4];
    res[5] = a[5];
    res[6] = a[6];
    res[7] = a[7];
    res[8] = a[8];
    res[9] = a[9];
    res[10] = a[10];
    res[11] = a[11];
    res[12] = a[12];
    res[13] = a[13];
    res[14] = a[14];
    res[15] = a[15];
    return res;
};

/// A rotor (Euclidean line) and translator (Ideal line)
template <typename ScalarType = float>
static PGA3D<ScalarType> rotor(float angle, PGA3D<ScalarType> line)
{
    return std::cos(angle / 2.0f) + std::sin(angle / 2.0f) * line.normalized();
}

/// Translator over directed line
template <typename ScalarType = float>
static PGA3D<ScalarType> translator(ScalarType dist, PGA3D<ScalarType> line)
{
    return 1.0f + dist / 2.0f * line;
}

/// A plane is defined using its homogenous equation ax + by + cz + d = 0
template <typename ScalarType = float>
static PGA3D<ScalarType> plane(ScalarType a, ScalarType b, ScalarType c, ScalarType d)
{
    /// PGA is plane based. Vectors are planes
    const PGA3D<ScalarType> e0(kE0), e1(kE1), e2(kE2), e3(kE3);
    return a * e1 + b * e2 + c * e3 + d * e0;
}

/// A point is just a homogeneous point, euclidean coordinates plus the origin
template <typename ScalarType = float>
static PGA3D<ScalarType> point(const ScalarType& x, const ScalarType& y, const ScalarType& z)
{
    const PGA3D<ScalarType> e0(kE0), e1(kE1), e2(kE2), e3(kE3);
    /// PGA points are trivectors.
    const PGA3D<ScalarType> e123 = e1 ^ e2 ^ e3, e032 = e0 ^ e3 ^ e2, e013 = e0 ^ e1 ^ e3, e021 = e0 ^ e2 ^ e1;
    return e123 + x * e032 + y * e013 + z * e021;
}

// for our toy problem (generate points on the surface of a torus)
// we start with a function that generates motors.
// circle(t) with t going from 0 to 1.
template <typename ScalarType = float>
static PGA3D<ScalarType> circle(ScalarType t, ScalarType radius, PGA3D<ScalarType> line)
{
    const PGA3D<ScalarType> e0(kE0), e1(kE1), e2(kE2), e3(kE3);
    return rotor(t * 2.0f * PI, line) * translator(radius, e1 * e0);
}

// a torus is now the product of two circles.
template <typename ScalarType = float>
static PGA3D<ScalarType> torus(ScalarType s,
                               ScalarType t,
                               ScalarType r1,
                               PGA3D<ScalarType> l1,
                               ScalarType r2,
                               PGA3D<ScalarType> l2)
{
    return circle(s, r2, l2) * circle(t, r1, l1);
}

// and to sample its points we simply sandwich the origin ..
template <typename ScalarType = float>
static PGA3D<ScalarType> point_on_torus(ScalarType s, ScalarType t)
{
    const PGA3D<ScalarType> e0(kE0), e1(kE1), e2(kE2), e3(kE3);
    const PGA3D<ScalarType> e123 = e1 ^ e2 ^ e3;
    PGA3D<ScalarType> to = torus(s, t, 0.25f, e1 * e2, 0.6f, e1 * e3);
    return to * e123 * ~to;
}

template <typename ScalarType = float>
void log(const PGA3D<ScalarType>& a, std::string name = "")
{
    static const char* basis[] = {
        "1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e12", "e31", "e23", "e021", "e013", "e032", "e123", "e0123"};

    if (!name.empty())
    {
        printf("%s = ", name.c_str());
    }
    int n = 0;
    for (int i = 0, j = 0; i < 16; i++)
        if (float(a[i]) != 0.0f)
        {
            n++;
            printf("%s%0.7g%s", (j > 0) ? " + " : "", float(a[i]), (i == 0) ? "" : basis[i]);
            j++;
        };
    if (n == 0)
        printf("0");
    printf("\n");
}

template <typename ScalarType = float>
PGA3D<ScalarType> motor_from_3_points_pairs(const std::array<PGA3D<ScalarType>, 3>& reference_points,
                                            const std::array<PGA3D<ScalarType>, 3>& target_points)
{
    const PGA3D<>& A = reference_points[0];
    const PGA3D<>& B = reference_points[1];
    const PGA3D<>& C = reference_points[2];
    const PGA3D<> A1 = target_points[0];
    const PGA3D<> B1 = target_points[1];
    const PGA3D<> C1 = target_points[2];

    const auto Va = (A1 * ~A).sqrt();
    const auto Ba = Va * B * ~Va;
    const auto Vb = ((A1 & B1) * ~(A1 & Ba)).sqrt();
    const auto Cba = Vb * Va * C * ~Va * ~Vb;
    const auto Vc = ((A1 & B1 & C1) * ~(A1 & B1 & Cba)).sqrt();

    return Vc * Vb * Va;
}

namespace float_basis
{
const static PGA3D<float> e0(kE0), e1(kE1), e2(kE2), e3(kE3);
const static PGA3D<float> e01(kE01), e02(kE02), e03(kE03);
const static PGA3D<float> e12(kE12), e23(kE23), e31(kE31);
const static PGA3D<float> e021(kE021), e023(kE032), e013(kE013), e123(kE123);
const static PGA3D<float> e0123(kE0123);
};