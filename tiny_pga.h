/*
 * This file is part of the Tiny-PGA distribution (https://github.com/sergehog/tiny_pga)
 * Copyright (c) 2020-2021 Sergey Smirnov / Seregium Oy.
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

//
// C++ implementation of 3D Projective Geometric Algebra (a.k.a Plane-based Geometric Algebra)
// C++ Templating helps with number of issues:
// * reducing computational complexity of PGA operators via compile-time optimizations
// * reducing memory footprint
// * compile-time checks for assignments correctness / blade matching (works to some extent)
//

#ifndef TINY_PGA_H
#define TINY_PGA_H

#include <array>
#include <cstdint>
#include <type_traits>

namespace tiny_pga
{

/// Bitmap of the presence of elements in the multivector
using Elems = std::uint16_t;

/// Internal namespace controlling compile-time operations with multivector elements
namespace elems
{

/// Elements are enumerated not in the logical order, but rather in the order they appear in the memory
enum class BitValues : Elems
{
    kScalar = (1U << 0U),
    // Vector
    kE0 = (1U << 1U),
    kE1 = (1U << 2U),
    kE2 = (1U << 3U),
    kE3 = (1U << 4U),
    // BivectorE
    kE12 = (1U << 5U),
    kE31 = (1U << 6U),
    kE23 = (1U << 7U),
    // Bivector0
    kE01 = (1U << 8U),
    kE02 = (1U << 9U),
    kE03 = (1U << 10U),
    // Trivector
    kE021 = (1U << 11U),
    kE013 = (1U << 12U),
    kE032 = (1U << 13U),
    kE123 = (1U << 14U),
    // Pseudo Scalar
    kE0123 = (1U << 15U),
};

constexpr bool has_scalar(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kScalar);
}

constexpr bool has_e0(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE0);
}

constexpr bool has_e1(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE1);
}

constexpr bool has_e2(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE2);
}

constexpr bool has_e3(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE3);
}

constexpr bool has_e01(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE01);
}

constexpr bool has_e02(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE02);
}

constexpr bool has_e03(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE03);
}

constexpr bool has_e12(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE12);
}

constexpr bool has_e31(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE31);
}

constexpr bool has_e23(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE23);
}

constexpr bool has_e021(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE021);
}

constexpr bool has_e013(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE013);
}

constexpr bool has_e032(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE032);
}

constexpr bool has_e123(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE123);
}

constexpr bool has_e0123(const Elems elem)
{
    return elem & static_cast<Elems>(BitValues::kE0123);
}

//constexpr bool has_vector(Elems elems)
//{
//    return has_e0(elems) || has_e1(elems) || has_e2(elems) || has_e3(elems);
//}
//
//constexpr bool has_bivectorE(Elems elems)
//{
//    return has_e23(elems) || has_e31(elems) || has_e12(elems) || has_scalar(elems);
//}
//
//constexpr bool has_bivector0(Elems elems)
//{
//    return has_e01(elems) || has_e02(elems) || has_e03(elems) || has_e0123(elems);
//}
//
//constexpr bool has_trivector(Elems elems)
//{
//    return has_e021(elems) || has_e013(elems) || has_e032(elems) || has_e123(elems);
//}

constexpr Elems elements(const bool scalar,
                         const bool e0,
                         const bool e1,
                         const bool e2,
                         const bool e3,
                         const bool e01,
                         const bool e02,
                         const bool e03,
                         const bool e12,
                         const bool e31,
                         const bool e23,
                         const bool e021,
                         const bool e013,
                         const bool e032,
                         const bool e123,
                         const bool e0123)
{
    Elems out_elements{};
    out_elements |= scalar ? static_cast<Elems>(BitValues::kScalar) : 0U;
    out_elements |= e0 ? static_cast<Elems>(BitValues::kE0) : 0U;
    out_elements |= e1 ? static_cast<Elems>(BitValues::kE1) : 0U;
    out_elements |= e2 ? static_cast<Elems>(BitValues::kE2) : 0U;
    out_elements |= e3 ? static_cast<Elems>(BitValues::kE3) : 0U;
    out_elements |= e01 ? static_cast<Elems>(BitValues::kE01) : 0U;
    out_elements |= e02 ? static_cast<Elems>(BitValues::kE02) : 0U;
    out_elements |= e03 ? static_cast<Elems>(BitValues::kE03) : 0U;
    out_elements |= e12 ? static_cast<Elems>(BitValues::kE12) : 0U;
    out_elements |= e31 ? static_cast<Elems>(BitValues::kE31) : 0U;
    out_elements |= e23 ? static_cast<Elems>(BitValues::kE23) : 0U;
    out_elements |= e021 ? static_cast<Elems>(BitValues::kE021) : 0U;
    out_elements |= e013 ? static_cast<Elems>(BitValues::kE013) : 0U;
    out_elements |= e032 ? static_cast<Elems>(BitValues::kE032) : 0U;
    out_elements |= e123 ? static_cast<Elems>(BitValues::kE123) : 0U;
    out_elements |= e0123 ? static_cast<Elems>(BitValues::kE0123) : 0U;
    return out_elements;
}

/// compile-time multivector array size
constexpr size_t count(Elems elems)
{
    size_t result{};
    result = static_cast<size_t>(has_scalar(elems));
    result += static_cast<size_t>(has_e0(elems));
    result += static_cast<size_t>(has_e1(elems));
    result += static_cast<size_t>(has_e2(elems));
    result += static_cast<size_t>(has_e3(elems));
    result += static_cast<size_t>(has_e01(elems));
    result += static_cast<size_t>(has_e02(elems));
    result += static_cast<size_t>(has_e03(elems));
    result += static_cast<size_t>(has_e12(elems));
    result += static_cast<size_t>(has_e31(elems));
    result += static_cast<size_t>(has_e23(elems));
    result += static_cast<size_t>(has_e021(elems));
    result += static_cast<size_t>(has_e013(elems));
    result += static_cast<size_t>(has_e032(elems));
    result += static_cast<size_t>(has_e123(elems));
    result += static_cast<size_t>(has_e0123(elems));
    return result;
}

/// indexing within dynamic (compile-time) multivector array
constexpr size_t index_scalar(Elems elems)
{
    return 0UL;
}

constexpr size_t index_e0(Elems elems)
{
    return has_scalar(elems);
}

constexpr size_t index_e1(Elems elems)
{
    return index_e0(elems) + static_cast<size_t>(has_e0(elems));
}

constexpr size_t index_e2(Elems elems)
{
    return index_e1(elems) + static_cast<size_t>(has_e1(elems));
}

constexpr size_t index_e3(Elems elems)
{
    return index_e2(elems) + static_cast<size_t>(has_e2(elems));
}

constexpr size_t index_e01(Elems elems)
{
    return index_e3(elems) + static_cast<size_t>(has_e3(elems));
}

constexpr size_t index_e02(Elems elems)
{
    return index_e01(elems) + static_cast<size_t>(has_e01(elems));
}

constexpr size_t index_e03(Elems elems)
{
    return index_e02(elems) + static_cast<size_t>(has_e02(elems));
}

constexpr size_t index_e12(Elems elems)
{
    return index_e03(elems) + static_cast<size_t>(has_e03(elems));
}

constexpr size_t index_e31(Elems elems)
{
    return index_e12(elems) + static_cast<size_t>(has_e12(elems));
}

constexpr size_t index_e23(Elems elems)
{
    return index_e31(elems) + static_cast<size_t>(has_e31(elems));
}

constexpr size_t index_e021(Elems elems)
{
    return index_e23(elems) + static_cast<size_t>(has_e23(elems));
}

constexpr size_t index_e013(Elems elems)
{
    return index_e021(elems) + static_cast<size_t>(has_e021(elems));
}

constexpr size_t index_e032(Elems elems)
{
    return index_e013(elems) + static_cast<size_t>(has_e013(elems));
}

constexpr size_t index_e123(Elems elems)
{
    return index_e032(elems) + static_cast<size_t>(has_e032(elems));
}

constexpr size_t index_e0123(Elems elems)
{
    return index_e123(elems) + static_cast<size_t>(has_e123(elems));
}

constexpr Elems geometric_product(const Elems elems1, const Elems elems2)
{
    const bool scalar = (has_scalar(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e1(elems2)) ||
                        (has_e2(elems1) && has_e2(elems2)) || (has_e3(elems1) && has_e3(elems2)) ||
                        (has_e12(elems1) && has_e12(elems2)) || (has_e31(elems1) && has_e31(elems2)) ||
                        (has_e23(elems1) && has_e23(elems2)) || (has_e123(elems1) && has_e123(elems2));

    const bool e0 = (has_scalar(elems1) && has_e0(elems2)) || (has_e0(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e1(elems2)) ||
                    (has_e2(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e2(elems2)) ||
                    (has_e3(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e3(elems2)) ||
                    (has_e12(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e12(elems2)) ||
                    (has_e31(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e31(elems2)) ||
                    (has_e23(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e23(elems2)) ||
                    (has_e123(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e123(elems2));

    const bool e1 = (has_scalar(elems1) && has_e1(elems2)) || (has_e1(elems1) && has_scalar(elems2)) ||
                    (has_e2(elems1) && has_e12(elems2)) || (has_e3(elems1) && has_e31(elems2)) ||
                    (has_e12(elems1) && has_e2(elems2)) || (has_e31(elems1) && has_e3(elems2)) ||
                    (has_e23(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e23(elems2));

    const bool e2 = (has_scalar(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e12(elems2)) || (has_e3(elems1) && has_e23(elems2)) ||
                    (has_e12(elems1) && has_e1(elems2)) || (has_e23(elems1) && has_e3(elems2)) ||
                    (has_e31(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e31(elems2));

    const bool e3 = (has_scalar(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e31(elems2)) || (has_e2(elems1) && has_e23(elems2)) ||
                    (has_e31(elems1) && has_e1(elems2)) || (has_e23(elems1) && has_e2(elems2)) ||
                    (has_e12(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e12(elems2));

    const bool e01 = (has_scalar(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e1(elems2)) || (has_e1(elems1) && has_e0(elems2)) ||
                     (has_e2(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e2(elems2)) ||
                     (has_e3(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e3(elems2)) ||
                     (has_e02(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e02(elems2)) ||
                     (has_e03(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e03(elems2)) ||
                     (has_e23(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e23(elems2)) ||
                     (has_e032(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e032(elems2));

    const bool e02 = (has_scalar(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_e0(elems2)) ||
                     (has_e1(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e1(elems2)) ||
                     (has_e3(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e3(elems2)) ||
                     (has_e01(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e01(elems2)) ||
                     (has_e03(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e03(elems2)) ||
                     (has_e31(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e31(elems2)) ||
                     (has_e013(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e013(elems2));

    const bool e03 = (has_scalar(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e0(elems2)) ||
                     (has_e1(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e1(elems2)) ||
                     (has_e2(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e2(elems2)) ||
                     (has_e01(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e01(elems2)) ||
                     (has_e02(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e02(elems2)) ||
                     (has_e12(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e12(elems2)) ||
                     (has_e021(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e021(elems2));

    const bool e12 = (has_scalar(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_e1(elems2)) ||
                     (has_e3(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e3(elems2)) ||
                     (has_e31(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e31(elems2));

    const bool e31 = (has_scalar(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e1(elems2)) ||
                     (has_e2(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e2(elems2)) ||
                     (has_e12(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e12(elems2));

    const bool e23 = (has_scalar(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_scalar(elems2)) ||
                     (has_e2(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e2(elems2)) ||
                     (has_e1(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e1(elems2)) ||
                     (has_e12(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e12(elems2));

    const bool e021 = (has_scalar(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e0(elems2)) ||
                      (has_e1(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e1(elems2)) ||
                      (has_e2(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e2(elems2)) ||
                      (has_e3(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e3(elems2)) ||
                      (has_e03(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e03(elems2)) ||
                      (has_e31(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e31(elems2)) ||
                      (has_e23(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e23(elems2));

    const bool e013 = (has_scalar(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e0(elems2)) ||
                      (has_e1(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e1(elems2)) ||
                      (has_e3(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e3(elems2)) ||
                      (has_e2(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e2(elems2)) ||
                      (has_e02(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e02(elems2)) ||
                      (has_e12(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e12(elems2)) ||
                      (has_e23(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e23(elems2));

    const bool e032 = (has_scalar(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e0(elems2)) ||
                      (has_e2(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e2(elems2)) ||
                      (has_e3(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e3(elems2)) ||
                      (has_e1(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e1(elems2)) ||
                      (has_e01(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e01(elems2)) ||
                      (has_e12(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e12(elems2)) ||
                      (has_e31(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e31(elems2));

    const bool e123 = (has_scalar(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_scalar(elems2)) ||
                      (has_e1(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e1(elems2)) ||
                      (has_e2(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e2(elems2)) ||
                      (has_e3(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e3(elems2));

    const bool e0123 = (has_scalar(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_scalar(elems2)) ||
                       (has_e0(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e0(elems2)) ||
                       (has_e1(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e1(elems2)) ||
                       (has_e2(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e2(elems2)) ||
                       (has_e3(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e3(elems2)) ||
                       (has_e01(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e01(elems2)) ||
                       (has_e02(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e02(elems2)) ||
                       (has_e03(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e03(elems2));

    return elements(scalar, e0, e1, e2, e3, e01, e02, e03, e12, e31, e23, e021, e013, e032, e123, e0123);
}

constexpr Elems addition(const Elems elems1, const Elems elems2)
{
    const bool scalar = has_scalar(elems1) || has_scalar(elems2);
    const bool e0 = has_e0(elems1) || has_e0(elems2);
    const bool e1 = has_e1(elems1) || has_e1(elems2);
    const bool e2 = has_e2(elems1) || has_e2(elems2);
    const bool e3 = has_e3(elems1) || has_e3(elems2);
    const bool e01 = has_e01(elems1) || has_e01(elems2);
    const bool e02 = has_e02(elems1) || has_e02(elems2);
    const bool e03 = has_e03(elems1) || has_e03(elems2);
    const bool e12 = has_e12(elems1) || has_e12(elems2);
    const bool e31 = has_e31(elems1) || has_e31(elems2);
    const bool e23 = has_e23(elems1) || has_e23(elems2);
    const bool e021 = has_e021(elems1) || has_e021(elems2);
    const bool e013 = has_e013(elems1) || has_e013(elems2);
    const bool e032 = has_e032(elems1) || has_e032(elems2);
    const bool e123 = has_e123(elems1) || has_e123(elems2);
    const bool e0123 = has_e0123(elems1) || has_e0123(elems2);

    return elements(scalar, e0, e1, e2, e3, e01, e02, e03, e12, e31, e23, e021, e013, e032, e123, e0123);
}

constexpr Elems inverse(Elems elems)
{
    Elems result{0U};
    result += has_e0123(elems) ? static_cast<Elems>(elems::BitValues::kScalar) : 0UL;
    result += has_e123(elems) ? static_cast<Elems>(elems::BitValues::kE0) : 0UL;
    result += has_e032(elems) ? static_cast<Elems>(elems::BitValues::kE1) : 0UL;
    result += has_e013(elems) ? static_cast<Elems>(elems::BitValues::kE2) : 0UL;
    result += has_e021(elems) ? static_cast<Elems>(elems::BitValues::kE3) : 0UL;

    result += has_e23(elems) ? static_cast<Elems>(elems::BitValues::kE01) : 0UL;
    result += has_e31(elems) ? static_cast<Elems>(elems::BitValues::kE02) : 0UL;
    result += has_e12(elems) ? static_cast<Elems>(elems::BitValues::kE03) : 0UL;

    result += has_e01(elems) ? static_cast<Elems>(elems::BitValues::kE23) : 0UL;
    result += has_e02(elems) ? static_cast<Elems>(elems::BitValues::kE31) : 0UL;
    result += has_e03(elems) ? static_cast<Elems>(elems::BitValues::kE12) : 0UL;

    result += has_e0(elems) ? static_cast<Elems>(elems::BitValues::kE123) : 0UL;
    result += has_e1(elems) ? static_cast<Elems>(elems::BitValues::kE032) : 0UL;
    result += has_e2(elems) ? static_cast<Elems>(elems::BitValues::kE013) : 0UL;
    result += has_e3(elems) ? static_cast<Elems>(elems::BitValues::kE021) : 0UL;

    result += has_scalar(elems) ? static_cast<Elems>(elems::BitValues::kE0123) : 0UL;
    return result;
}

constexpr Elems PlaneElems = static_cast<Elems>(elems::BitValues::kE0) + static_cast<Elems>(elems::BitValues::kE1) +
                             static_cast<Elems>(elems::BitValues::kE2) + static_cast<Elems>(elems::BitValues::kE3);

constexpr Elems LineElems = static_cast<Elems>(elems::BitValues::kE01) + static_cast<Elems>(elems::BitValues::kE02) +
                            static_cast<Elems>(elems::BitValues::kE03) + static_cast<Elems>(elems::BitValues::kE23) +
                            static_cast<Elems>(elems::BitValues::kE31) + static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems PointElems = static_cast<Elems>(elems::BitValues::kE123) + static_cast<Elems>(elems::BitValues::kE032) +
                             static_cast<Elems>(elems::BitValues::kE013) + static_cast<Elems>(elems::BitValues::kE021);

constexpr Elems RotorElems = static_cast<Elems>(elems::BitValues::kScalar) +
                             static_cast<Elems>(elems::BitValues::kE23) + static_cast<Elems>(elems::BitValues::kE31) +
                             static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems TranslatorElems =
    static_cast<Elems>(elems::BitValues::kE01) + static_cast<Elems>(elems::BitValues::kE02) +
    static_cast<Elems>(elems::BitValues::kE03) + static_cast<Elems>(elems::BitValues::kScalar);

constexpr Elems MotorElems = static_cast<Elems>(elems::BitValues::kScalar) +
                             static_cast<Elems>(elems::BitValues::kE23) + static_cast<Elems>(elems::BitValues::kE31) +
                             static_cast<Elems>(elems::BitValues::kE12) + static_cast<Elems>(elems::BitValues::kE01) +
                             static_cast<Elems>(elems::BitValues::kE02) + static_cast<Elems>(elems::BitValues::kE03) +
                             static_cast<Elems>(elems::BitValues::kE0123);

}  // namespace elems

template <Elems elements, typename ScalarType = float>
struct Multivector;
using PlaneF = Multivector<elems::PlaneElems, float>;
using LineF = Multivector<elems::LineElems, float>;
using PointF = Multivector<elems::PointElems, float>;
using RotorF = Multivector<elems::RotorElems, float>;
using TranslatorF = Multivector<elems::TranslatorElems, float>;
using MotorF = Multivector<elems::MotorElems, float>;

/// Compile-time optimized (using templating) implementation of 3D PGA Multivector
template <Elems elements, typename ScalarType>
struct Multivector
{
    using Plane = Multivector<elems::PlaneElems, ScalarType>;
    using Line = Multivector<elems::LineElems, ScalarType>;
    using Point = Multivector<elems::PointElems, ScalarType>;
    using Rotor = Multivector<elems::RotorElems, ScalarType>;
    using Translator = Multivector<elems::TranslatorElems, ScalarType>;
    using Motor = Multivector<elems::MotorElems, ScalarType>;
    /// Exposing elements outside, as static const value
    static constexpr Elems Elements = elements;

    /// Actual multivector values
    std::array<ScalarType, elems::count(elements)> values;

#define ELEM_ACCESS_FUNCTION(element_name)                                                      \
    template <typename T = ScalarType>                                                          \
    typename std::enable_if<elems::has_##element_name(elements), T>::type& element_name()       \
    {                                                                                           \
        return values[elems::index_##element_name(elements)];                                   \
    };                                                                                          \
    template <typename T = ScalarType>                                                          \
    typename std::enable_if<elems::has_##element_name(elements), T>::type element_name() const  \
    {                                                                                           \
        return values[elems::index_##element_name(elements)];                                   \
    }                                                                                           \
    template <typename T = ScalarType>                                                          \
    typename std::enable_if<!elems::has_##element_name(elements), T>::type& element_name()      \
    {                                                                                           \
        static ScalarType stub_element{};                                                       \
        stub_element = 0;                                                                       \
        return stub_element;                                                                    \
    }                                                                                           \
    template <typename T = ScalarType>                                                          \
    typename std::enable_if<!elems::has_##element_name(elements), T>::type element_name() const \
    {                                                                                           \
        return 0.;                                                                              \
    }

    ELEM_ACCESS_FUNCTION(scalar);

    ELEM_ACCESS_FUNCTION(e0);
    ELEM_ACCESS_FUNCTION(e1);
    ELEM_ACCESS_FUNCTION(e2);
    ELEM_ACCESS_FUNCTION(e3);

    ELEM_ACCESS_FUNCTION(e01);
    ELEM_ACCESS_FUNCTION(e02);
    ELEM_ACCESS_FUNCTION(e03);

    ELEM_ACCESS_FUNCTION(e12);
    ELEM_ACCESS_FUNCTION(e31);
    ELEM_ACCESS_FUNCTION(e23);

    ELEM_ACCESS_FUNCTION(e021);
    ELEM_ACCESS_FUNCTION(e013);
    ELEM_ACCESS_FUNCTION(e032);

    ELEM_ACCESS_FUNCTION(e123);
    ELEM_ACCESS_FUNCTION(e0123);

    /// Geometric Product operator
    /// Resulting multivector type deducted automatically using templating
    template <Elems other_elements, typename ScalarTypeIn>
    Multivector<elems::geometric_product(elements, other_elements), ScalarType> operator*(
        const Multivector<other_elements, ScalarTypeIn>& other) const
    {
        constexpr Elems out_elems = elems::geometric_product(elements, other_elements);
        Multivector<out_elems, ScalarType> out{};

#define ELEM_MULTIPLY(elem_out, elem1, elem2, sign1)                        \
    if (elems::has_##elem1(elements) && elems::has_##elem2(other_elements)) \
    {                                                                       \
        out.elem_out() sign1 elem1() * other.elem2();                       \
    }

#define ELEM_BOTH_MULTIPLY(elem_out, elem1, elem2, sign1, sign2)            \
    if (elems::has_##elem1(elements) && elems::has_##elem2(other_elements)) \
    {                                                                       \
        out.elem_out() sign1 elem1() * other.elem2();                       \
    }                                                                       \
    if (elems::has_##elem2(elements) && elems::has_##elem1(other_elements)) \
    {                                                                       \
        out.elem_out() sign2 elem2() * other.elem1();                       \
    }

        if (elems::has_scalar(out_elems))
        {
            ELEM_MULTIPLY(scalar, scalar, scalar, +=);
            ELEM_MULTIPLY(scalar, e1, e1, +=);
            ELEM_MULTIPLY(scalar, e2, e2, +=);
            ELEM_MULTIPLY(scalar, e3, e3, +=);
            ELEM_MULTIPLY(scalar, e12, e12, -=);
            ELEM_MULTIPLY(scalar, e31, e31, -=);
            ELEM_MULTIPLY(scalar, e23, e23, -=);
            ELEM_MULTIPLY(scalar, e123, e123, -=);
        }

        if (elems::has_e0(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e0, scalar, e0, +=, +=);
            ELEM_BOTH_MULTIPLY(e0, e1, e01, -=, +=);
            ELEM_BOTH_MULTIPLY(e0, e2, e02, -=, +=);
            ELEM_BOTH_MULTIPLY(e0, e3, e03, -=, +=);
            ELEM_BOTH_MULTIPLY(e0, e12, e021, +=, +=);
            ELEM_BOTH_MULTIPLY(e0, e31, e013, +=, +=);
            ELEM_BOTH_MULTIPLY(e0, e23, e032, +=, +=);
            ELEM_BOTH_MULTIPLY(e0, e123, e0123, +=, -=);
        }

        if (elems::has_e1(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e1, scalar, e1, +=, +=);
            ELEM_BOTH_MULTIPLY(e1, e2, e12, -=, +=);
            ELEM_BOTH_MULTIPLY(e1, e3, e31, +=, -=);
            ELEM_BOTH_MULTIPLY(e1, e23, e123, -=, -=);
        }

        if (elems::has_e2(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e2, scalar, e2, +=, +=);
            ELEM_BOTH_MULTIPLY(e2, e1, e12, +=, -=);
            ELEM_BOTH_MULTIPLY(e2, e3, e23, -=, +=);
            ELEM_BOTH_MULTIPLY(e2, e31, e123, -=, -=);
        }

        if (elems::has_e3(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e3, scalar, e3, +=, +=);
            ELEM_BOTH_MULTIPLY(e3, e1, e31, -=, +=);
            ELEM_BOTH_MULTIPLY(e3, e2, e23, +=, -=);
            ELEM_BOTH_MULTIPLY(e3, e12, e123, -=, -=);
        }

        if (elems::has_e01(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e01, scalar, e01, +=, +=);
            ELEM_BOTH_MULTIPLY(e01, e0, e1, +=, -=);
            ELEM_BOTH_MULTIPLY(e01, e2, e021, -=, -=);
            ELEM_BOTH_MULTIPLY(e01, e3, e013, +=, +=);
            ELEM_BOTH_MULTIPLY(e01, e02, e12, -=, +=);
            ELEM_BOTH_MULTIPLY(e01, e03, e31, +=, -=);
            ELEM_BOTH_MULTIPLY(e01, e23, e0123, -=, -=);
            ELEM_BOTH_MULTIPLY(e01, e032, e123, +=, -=);
        }

        if (elems::has_e02(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e02, scalar, e02, +=, +=);
            ELEM_BOTH_MULTIPLY(e02, e0, e2, +=, -=);
            ELEM_BOTH_MULTIPLY(e02, e1, e021, +=, +=);
            ELEM_BOTH_MULTIPLY(e02, e3, e032, -=, -=);
            ELEM_BOTH_MULTIPLY(e02, e01, e12, +=, -=);
            ELEM_BOTH_MULTIPLY(e02, e03, e23, -=, +=);
            ELEM_BOTH_MULTIPLY(e02, e31, e0123, -=, -=);
            ELEM_BOTH_MULTIPLY(e02, e013, e123, +=, -=);
        }

        if (elems::has_e03(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e03, scalar, e03, +=, +=);
            ELEM_BOTH_MULTIPLY(e03, e0, e3, +=, -=);
            ELEM_BOTH_MULTIPLY(e03, e1, e013, -=, -=);
            ELEM_BOTH_MULTIPLY(e03, e2, e032, +=, +=);
            ELEM_BOTH_MULTIPLY(e03, e01, e31, -=, +=);
            ELEM_BOTH_MULTIPLY(e03, e02, e23, +=, -=);
            ELEM_BOTH_MULTIPLY(e03, e12, e0123, -=, -=);
            ELEM_BOTH_MULTIPLY(e03, e021, e123, +=, -=);
        }

        if (elems::has_e12(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e12, scalar, e12, +=, +=);
            ELEM_BOTH_MULTIPLY(e12, e1, e2, +=, -=);
            ELEM_BOTH_MULTIPLY(e12, e3, e123, +=, +=);
            ELEM_BOTH_MULTIPLY(e12, e31, e23, +=, -=);
        }

        if (elems::has_e31(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e31, scalar, e31, +=, +=);
            ELEM_BOTH_MULTIPLY(e31, e3, e1, +=, -=);
            ELEM_BOTH_MULTIPLY(e31, e2, e123, +=, +=);
            ELEM_BOTH_MULTIPLY(e31, e12, e23, -=, +=);
        }

        if (elems::has_e23(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e23, scalar, e23, +=, +=);
            ELEM_BOTH_MULTIPLY(e23, e2, e3, +=, -=);
            ELEM_BOTH_MULTIPLY(e23, e1, e123, +=, +=);
            ELEM_BOTH_MULTIPLY(e23, e12, e31, +=, -=);
        }

        if (elems::has_e021(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e021, scalar, e021, +=, +=);
            ELEM_BOTH_MULTIPLY(e021, e0, e12, +=, -=);
            ELEM_BOTH_MULTIPLY(e021, e1, e02, +=, +=);
            ELEM_BOTH_MULTIPLY(e021, e2, e01, -=, -=);
            ELEM_BOTH_MULTIPLY(e021, e3, e0123, +=, -=);
            ELEM_BOTH_MULTIPLY(e021, e03, e123, +=, -=);
            ELEM_BOTH_MULTIPLY(e021, e31, e032, +=, -=);
            ELEM_BOTH_MULTIPLY(e021, e23, e013, -=, +=);
        }

        if (elems::has_e013(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e013, scalar, e013, +=, +=);
            ELEM_BOTH_MULTIPLY(e013, e0, e31, +=, -=);
            ELEM_BOTH_MULTIPLY(e013, e1, e03, +=, -=);
            ELEM_BOTH_MULTIPLY(e013, e3, e01, +=, +=);
            ELEM_BOTH_MULTIPLY(e013, e2, e0123, +=, -=);
            ELEM_BOTH_MULTIPLY(e013, e02, e123, +=, -=);
            ELEM_BOTH_MULTIPLY(e013, e12, e032, -=, +=);
            ELEM_BOTH_MULTIPLY(e013, e23, e021, -=, +=);
        }

        if (elems::has_e032(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e032, scalar, e032, +=, +=);
            ELEM_BOTH_MULTIPLY(e032, e0, e23, -=, -=);
            ELEM_BOTH_MULTIPLY(e032, e2, e03, +=, +=);
            ELEM_BOTH_MULTIPLY(e032, e3, e02, -=, -=);
            ELEM_BOTH_MULTIPLY(e032, e1, e0123, +=, -=);
            ELEM_BOTH_MULTIPLY(e032, e01, e123, -=, +=);
            ELEM_BOTH_MULTIPLY(e032, e12, e013, +=, -=);
            ELEM_BOTH_MULTIPLY(e032, e31, e021, -=, +=);
        }

        if (elems::has_e123(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e123, scalar, e123, +=, +=);
            ELEM_BOTH_MULTIPLY(e123, e1, e23, +=, +=);
            ELEM_BOTH_MULTIPLY(e123, e2, e31, +=, +=);
            ELEM_BOTH_MULTIPLY(e123, e3, e12, +=, +=);
        }

        if (elems::has_e0123(out_elems))
        {
            ELEM_BOTH_MULTIPLY(e0123, scalar, e0123, +=, +=);
            ELEM_BOTH_MULTIPLY(e0123, e0, e123, +=, -=);
            ELEM_BOTH_MULTIPLY(e0123, e1, e032, +=, -=);
            ELEM_BOTH_MULTIPLY(e0123, e2, e013, +=, -=);
            ELEM_BOTH_MULTIPLY(e0123, e3, e021, +=, -=);
            ELEM_BOTH_MULTIPLY(e0123, e01, e23, +=, +=);
            ELEM_BOTH_MULTIPLY(e0123, e02, e31, +=, +=);
            ELEM_BOTH_MULTIPLY(e0123, e03, e12, +=, +=);
        }
#undef ELEM_BOTH_MULTIPLY
        return out;
    }

    /// Multivector addition
    template <Elems other_elements, typename ScalarTypeIn>
    Multivector<elems::addition(elements, other_elements), ScalarType> operator+(
        const Multivector<other_elements, ScalarTypeIn>& other) const
    {
        constexpr Elems out_elems = elems::addition(elements, other_elements);
        Multivector<out_elems, ScalarType> out{};
#define ELEM_ADD(name) \
        if(elems::has_##name(out_elems)) \
        { \
            if(elems::has_##name(elements)) \
            { \
                out.name() = name(); \
            } \
            if(elems::has_##name(other_elements)) \
            { \
                out.name() += other.name(); \
            } \
        }

        ELEM_ADD(scalar);
        ELEM_ADD(e0);
        ELEM_ADD(e1);
        ELEM_ADD(e2);
        ELEM_ADD(e3);
        ELEM_ADD(e01);
        ELEM_ADD(e02);
        ELEM_ADD(e03);
        ELEM_ADD(e12);
        ELEM_ADD(e31);
        ELEM_ADD(e23);
        ELEM_ADD(e021);
        ELEM_ADD(e013);
        ELEM_ADD(e032);
        ELEM_ADD(e123);
        ELEM_ADD(e0123);
        return out;
    }

    Multivector<elements, ScalarType> reverse() const
    {
        Multivector<elements, ScalarType> out{};
        if(elems::has_scalar(elements))
        {
            out.scalar() = scalar();
        }
        if(elems::has_e0(elements))
        {
            out.e0() = e0();
        }
        if(elems::has_e1(elements))
        {
            out.e1() = e1();
        }
        if(elems::has_e2(elements))
        {
            out.e2() = e2();
        }
        if(elems::has_e3(elements))
        {
            out.e3() = e3();
        }
        if(elems::has_e01(elements))
        {
            out.e01() = -e01();
        }
        if(elems::has_e02(elements))
        {
            out.e02() = -e02();
        }
        if(elems::has_e03(elements))
        {
            out.e03() = -e03();
        }
        if(elems::has_e12(elements))
        {
            out.e12() = -e12();
        }
        if(elems::has_e31(elements))
        {
            out.e31() = -e31();
        }
        if(elems::has_e23(elements))
        {
            out.e23() = -e23();
        }
        if(elems::has_e021(elements))
        {
            out.e021() = -e021();
        }
        if(elems::has_e013(elements))
        {
            out.e013() = -e013();
        }
        if(elems::has_e032(elements))
        {
            out.e032() = -e032();
        }
        if(elems::has_e123(elements))
        {
            out.e123() = -e123();
        }
        if(elems::has_e0123(elements))
        {
            out.e1e012323() = e0123();
        }
    }

    /// Reverse operator
    Multivector<elements, ScalarType> operator~() const
    {
        return reverse();
    }

    Multivector<elems::inverse(elements), ScalarType> inverse()
    {
        Multivector<elems::inverse(elements), ScalarType> result;
        if (elems::has_e0123(elements))
        {
            result.scalar() = e0123();
        }
        if (elems::has_e123(elements))
        {
            result.e0() = e123();
        }
        if (elems::has_e032(elements))
        {
            result.e1() = e032();
        }
        if (elems::has_e013(elements))
        {
            result.e2() = e013();
        }
        if (elems::has_e021(elements))
        {
            result.e3() = e021();
        }
        if (elems::has_e23(elements))
        {
            result.e01() = e23();
        }
        if (elems::has_e31(elements))
        {
            result.e02() = e31();
        }
        if (elems::has_e12(elements))
        {
            result.e03() = e12();
        }
        if (elems::has_e01(elements))
        {
            result.e23() = e01();
        }
        if (elems::has_e02(elements))
        {
            result.e31() = e02();
        }
        if (elems::has_e03(elements))
        {
            result.e12() = e03();
        }

        if (elems::has_e0(elements))
        {
            result.e123() = e0();
        }
        if (elems::has_e1(elements))
        {
            result.e032() = e1();
        }
        if (elems::has_e2(elements))
        {
            result.e013() = e2();
        }
        if (elems::has_e3(elements))
        {
            result.e021() = e3();
        }

        if (elems::has_scalar(elements))
        {
            result.e0123() = scalar();
        }
    }

    //    Multivector<elements, ScalarType> sandwich(const Motor& motor) const
    //    {
    //        Multivector<elements, ScalarType> out;
    //        auto res = motor * (*this) * (~motor);
    //        if (elems::has_vector(elements))
    //        {
    //            out.Vector.value = res.Vector.value;
    //        }
    //        if (elems::has_bivector0(elements))
    //        {
    //            out.Bivector0.value = res.Bivector0.value;
    //        }
    //        if (elems::has_bivectorE(elements))
    //        {
    //            out.BivectorE.value = res.BivectorE.value;
    //        }
    //        if (elems::has_trivector(elements))
    //        {
    //            out.Trivector.value = res.Trivector.value;
    //        }
    //        return out;
    //    }
    //
    //    Multivector<elements, ScalarType> sandwich(const Rotor& rotor) const
    //    {
    //        Multivector<elements, ScalarType> out;
    //        auto res = rotor * (*this) * (~rotor);
    //        if (elems::has_vector(elements))
    //        {
    //            out.Vector.value = res.Vector.value;
    //        }
    //        if (elems::has_bivector0(elements))
    //        {
    //            out.Bivector0.value = res.Bivector0.value;
    //        }
    //        if (elems::has_bivectorE(elements))
    //        {
    //            out.BivectorE.value = res.BivectorE.value;
    //        }
    //        if (elems::has_trivector(elements))
    //        {
    //            out.Trivector.value = res.Trivector.value;
    //        }
    //        return out;
    //    }
    //
    //    Multivector<elements, ScalarType> sandwich(const Translator& translator) const
    //    {
    //        Multivector<elements, ScalarType> out;
    //        auto res = translator * (*this) * (~translator);
    //        if (elems::has_vector(elements))
    //        {
    //            out.Vector.value = res.Vector.value;
    //        }
    //        if (elems::has_bivector0(elements))
    //        {
    //            out.Bivector0.value = res.Bivector0.value;
    //        }
    //        if (elems::has_bivectorE(elements))
    //        {
    //            out.BivectorE.value = res.BivectorE.value;
    //        }
    //        if (elems::has_trivector(elements))
    //        {
    //            out.Trivector.value = res.Trivector.value;
    //        }
    //        return out;
    //    }

    /// Converts whatever it's now to a plane
    explicit operator Plane() const
    {
        Multivector<elems::PlaneElems, ScalarType> plane{};
        if (elems::has_e0(elements))
        {
            plane.e0() = e0();
        }
        if (elems::has_e1(elements))
        {
            plane.e1() = e1();
        }
        if (elems::has_e2(elements))
        {
            plane.e2() = e2();
        }
        if (elems::has_e3(elements))
        {
            plane.e3() = e3();
        }
        return plane;
    }

    /// Converts whatever it's now to a point
    explicit operator Point() const
    {
        Multivector<elems::PointElems, ScalarType> point{};
        if (elems::has_e021(elements))
        {
            point.e021() = e021();
        }
        if (elems::has_e013(elements))
        {
            point.e013() = e013();
        }
        if (elems::has_e032(elements))
        {
            point.e032() = e032();
        }
        point.e123() = 1.F;
        return point;
    }
};
//
// template<template <Elems elements, typename ScalarType>
// struct Multivector<elements, ScalarType>
//{
//    Multivector()
//    {
//    }
//};
//
//
// template<>
// struct Multivector<elems::PointElems, float>
//{
//    Multivector(const float x, const float y, const float z)
//    {
//        Trivector.value = {x, y, z, 1.F};
//    }
//};

}  // namespace tiny_pga

#endif  // TINY_PGA_H
