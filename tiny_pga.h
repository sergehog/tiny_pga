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
// C++ Implementation of 3D Projectove Geometric Algebra (a.k.a Plane-based Geometric Algebra)
// Highly-templatized implementation helps with number of issues:
// * reducing computational complexity of PGA operators via compile-time optimizations
// * reducing memory footprint
// * compile-time checks for assignments correctness / blade matching (to some extent)
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

namespace elems
{

/// Elements are enumerated not in the logical order, but rather in the order they appear in the
/// memory
enum class BitValues : Elems
{
    // Vector
    kE0 = (1U << 0U),
    kE1 = (1U << 1U),
    kE2 = (1U << 2U),
    kE3 = (1U << 3U),
    // BivectorE
    kScalar = (1U << 4U),
    kE12 = (1U << 5U),
    kE31 = (1U << 6U),
    kE23 = (1U << 7U),
    // Bivector0
    kE01 = (1U << 8U),
    kE02 = (1U << 9U),
    kE03 = (1U << 10U),
    kE0123 = (1U << 11U),
    // Trivector
    kE021 = (1U << 12U),
    kE013 = (1U << 13U),
    kE032 = (1U << 14U),
    kE123 = (1U << 15U),
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

constexpr bool has_vector(Elems elems)
{
    return has_e0(elems) || has_e1(elems) || has_e2(elems) || has_e3(elems);
}

constexpr bool has_bivectorE(Elems elems)
{
    return has_e23(elems) || has_e31(elems) || has_e12(elems) || has_scalar(elems);
}

constexpr bool has_bivector0(Elems elems)
{
    return has_e01(elems) || has_e02(elems) || has_e03(elems) || has_e0123(elems);
}

constexpr bool has_trivector(Elems elems)
{
    return has_e021(elems) || has_e013(elems) || has_e032(elems) || has_e123(elems);
}

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
    out_elements |= scalar ? static_cast<Elems>(BitValues::kScalar) : 0;
    out_elements |= e0 ? static_cast<Elems>(BitValues::kE0) : 0;
    out_elements |= e1 ? static_cast<Elems>(BitValues::kE1) : 0;
    out_elements |= e2 ? static_cast<Elems>(BitValues::kE2) : 0;
    out_elements |= e3 ? static_cast<Elems>(BitValues::kE3) : 0;
    out_elements |= e01 ? static_cast<Elems>(BitValues::kE01) : 0;
    out_elements |= e02 ? static_cast<Elems>(BitValues::kE02) : 0;
    out_elements |= e03 ? static_cast<Elems>(BitValues::kE03) : 0;
    out_elements |= e12 ? static_cast<Elems>(BitValues::kE12) : 0;
    out_elements |= e31 ? static_cast<Elems>(BitValues::kE31) : 0;
    out_elements |= e23 ? static_cast<Elems>(BitValues::kE23) : 0;
    out_elements |= e021 ? static_cast<Elems>(BitValues::kE021) : 0;
    out_elements |= e013 ? static_cast<Elems>(BitValues::kE013) : 0;
    out_elements |= e032 ? static_cast<Elems>(BitValues::kE032) : 0;
    out_elements |= e123 ? static_cast<Elems>(BitValues::kE123) : 0;
    out_elements |= e0123 ? static_cast<Elems>(BitValues::kE0123) : 0;
    return out_elements;
}

constexpr Elems multiplication(const Elems elems1, const Elems elems2)
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

constexpr Elems PlaneElems = static_cast<Elems>(elems::BitValues::kE0) | static_cast<Elems>(elems::BitValues::kE1) |
                             static_cast<Elems>(elems::BitValues::kE2) | static_cast<Elems>(elems::BitValues::kE3);

constexpr Elems LineElems = static_cast<Elems>(elems::BitValues::kE01) | static_cast<Elems>(elems::BitValues::kE02) |
                            static_cast<Elems>(elems::BitValues::kE03) | static_cast<Elems>(elems::BitValues::kE23) |
                            static_cast<Elems>(elems::BitValues::kE31) | static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems PointElems = static_cast<Elems>(elems::BitValues::kE123) | static_cast<Elems>(elems::BitValues::kE032) |
                             static_cast<Elems>(elems::BitValues::kE013) | static_cast<Elems>(elems::BitValues::kE021);

constexpr Elems RotorElems = static_cast<Elems>(elems::BitValues::kScalar) |
                             static_cast<Elems>(elems::BitValues::kE23) | static_cast<Elems>(elems::BitValues::kE31) |
                             static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems TranslatorElems =
    static_cast<Elems>(elems::BitValues::kE01) | static_cast<Elems>(elems::BitValues::kE02) |
    static_cast<Elems>(elems::BitValues::kE03) | static_cast<Elems>(elems::BitValues::kScalar);

constexpr Elems MotorElems = static_cast<Elems>(elems::BitValues::kScalar) |
                             static_cast<Elems>(elems::BitValues::kE23) | static_cast<Elems>(elems::BitValues::kE31) |
                             static_cast<Elems>(elems::BitValues::kE12) | static_cast<Elems>(elems::BitValues::kE01) |
                             static_cast<Elems>(elems::BitValues::kE02) | static_cast<Elems>(elems::BitValues::kE03) |
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

    //    template <bool Condition, typename T>
    //    struct Conditional
    //    {
    //        T value;
    //    };
    //    template <typename T>
    //    struct Conditional<false, T>
    //    {
    //    };
    //
    //
    //    // Optimization of Memory Footprint with use of conditional elements
    //    Conditional<elems::has_vector(elements), std::array<ScalarType, 4U>> Vector;
    //    Conditional<elems::has_bivectorE(elements), std::array<ScalarType, 4U>> BivectorE;
    //    Conditional<elems::has_bivector0(elements), std::array<ScalarType, 4U>> Bivector0;
    //    Conditional<elems::has_trivector(elements), std::array<ScalarType, 4U>> Trivector;
    // In this macro we define setter and read-only getter functions,
    // as well as define private stub function, in case if element does not exist
    //#define ELEM_ACCESS_FUNCTION(element_name, array_position)                                      \
//    template <typename T = ScalarType>                                                          \
//    typename std::enable_if<elems::has_##element_name(elements), T>::type& element_name()       \
//    {                                                                                           \
//        return array_position;                                                                  \
//    };                                                                                          \
//    template <typename T = ScalarType>                                                          \
//    typename std::enable_if<elems::has_##element_name(elements), T>::type element_name() const  \
//    {                                                                                           \
//        return array_position;                                                                  \
//    }                                                                                           \
//    template <typename T = ScalarType>                                                          \
//    typename std::enable_if<!elems::has_##element_name(elements), T>::type& element_name()      \
//    {                                                                                           \
//        return stub_element;                                                                    \
//    }                                                                                           \
//    template <typename T = ScalarType>                                                          \
//    typename std::enable_if<!elems::has_##element_name(elements), T>::type element_name() const \
//    {                                                                                           \
//        return 0.;                                                                              \
//    }

    std::array<ScalarType, 4U> Vector;
    std::array<ScalarType, 4U> BivectorE;
    std::array<ScalarType, 4U> Bivector0;
    std::array<ScalarType, 4U> Trivector;

#define ELEM_ACCESS_FUNCTION(element_name, array_position) \
    ScalarType& element_name() { return array_position; }; \
    ScalarType element_name() const { return array_position; }

    ELEM_ACCESS_FUNCTION(e0, Vector[0]);
    ELEM_ACCESS_FUNCTION(e1, Vector[1]);
    ELEM_ACCESS_FUNCTION(e2, Vector[2]);
    ELEM_ACCESS_FUNCTION(e3, Vector[3]);

    ELEM_ACCESS_FUNCTION(e01, Bivector0[0]);
    ELEM_ACCESS_FUNCTION(e02, Bivector0[1]);
    ELEM_ACCESS_FUNCTION(e03, Bivector0[2]);
    ELEM_ACCESS_FUNCTION(e0123, Bivector0[3]);

    ELEM_ACCESS_FUNCTION(scalar, BivectorE[0]);
    ELEM_ACCESS_FUNCTION(e12, BivectorE[1]);
    ELEM_ACCESS_FUNCTION(e31, BivectorE[2]);
    ELEM_ACCESS_FUNCTION(e23, BivectorE[3]);

    ELEM_ACCESS_FUNCTION(e021, Trivector[0]);
    ELEM_ACCESS_FUNCTION(e013, Trivector[1]);
    ELEM_ACCESS_FUNCTION(e032, Trivector[2]);
    ELEM_ACCESS_FUNCTION(e123, Trivector[3]);
#undef DEFINE_ELEM_FUNCTION

    //    template <class T = Multivector<elems::RotorElems>>
    //    typename std::enable_if<elems::has_bivectorE(elements), T>::type rotor()
    //    {
    //        return Rotor{BivectorE.value};
    //    }
    //
    //    template <class T = Multivector<elems::TranslatorElems>>
    //    typename std::enable_if<elems::has_bivector0(elements), T>::type translator()
    //    {
    //        return Translator{Bivector0.value};
    //    }
    //
    //    template <class T = Multivector<elems::MotorElems>>
    //    typename std::enable_if<elems::has_bivector0(elements), T>::type motor()
    //    {
    //        return Motor{BivectorE.value, Bivector0.value};
    //    }

    /// Generic multiplication operator
    /// When using it, elements of resulting multivector (i.e. Blades) might grow, even though their
    /// real values remain zeros Consider casting result back to desired type in order to keep
    /// compile-time constraints active
    template <Elems other_elements, typename ScalarTypeIn>
    Multivector<elems::multiplication(elements, other_elements), ScalarType> operator*(
        const Multivector<other_elements, ScalarTypeIn>& other) const
    {
        constexpr Elems out_elems = elems::multiplication(elements, other_elements);
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

    /// Reverse operator
    Multivector<elements, ScalarType> operator~() const
    {
        Multivector<elements, ScalarType> out{};
        if (elems::has_vector(elements))
        {
            out.Vector[0] = Vector[0];
            out.Vector[1] = Vector[1];
            out.Vector[2] = Vector[2];
            out.Vector[3] = Vector[3];
        }

        if (elems::has_bivector0(elements))
        {
            out.Bivector0[0] = -Bivector0[0];
            out.Bivector0[1] = -Bivector0[1];
            out.Bivector0[2] = -Bivector0[2];
            out.Bivector0[3] = Bivector0[3];
        }

        if (elems::has_bivectorE(elements))
        {
            out.BivectorE[0] = BivectorE[0];
            out.BivectorE[1] = -BivectorE[1];
            out.BivectorE[2] = -BivectorE[2];
            out.BivectorE[3] = -BivectorE[3];
        }

        if (elems::has_trivector(elements))
        {
            out.Trivector[0] = -Trivector[0];
            out.Trivector[1] = -Trivector[1];
            out.Trivector[2] = -Trivector[2];
            out.Trivector[3] = -Trivector[3];
        }
        return out;
    }

    Multivector<elements, ScalarType> sandwich(const Motor& motor) const
    {
        Multivector<elements, ScalarType> out;
        auto res = motor * (*this) * (~motor);
        if (elems::has_vector(elements))
        {
            out.Vector = res.Vector;
        }
        if (elems::has_bivector0(elements))
        {
            out.Bivector0 = res.Bivector0;
        }
        if (elems::has_bivectorE(elements))
        {
            out.BivectorE = res.BivectorE;
        }
        if (elems::has_trivector(elements))
        {
            out.Trivector = res.Trivector;
        }
        return out;
    }

    Multivector<elements, ScalarType> sandwich(const Rotor& rotor) const
    {
        Multivector<elements, ScalarType> out;
        auto res = rotor * (*this) * (~rotor);
        if (elems::has_vector(elements))
        {
            out.Vector = res.Vector;
        }
        if (elems::has_bivector0(elements))
        {
            out.Bivector0 = res.Bivector0;
        }
        if (elems::has_bivectorE(elements))
        {
            out.BivectorE = res.BivectorE;
        }
        if (elems::has_trivector(elements))
        {
            out.Trivector = res.Trivector;
        }
        return out;
    }

    Multivector<elements, ScalarType> sandwich(const Translator& translator) const
    {
        Multivector<elements, ScalarType> out;
        auto res = translator * (*this) * (~translator);
        if (elems::has_vector(elements))
        {
            out.Vector = res.Vector;
        }
        if (elems::has_bivector0(elements))
        {
            out.Bivector0 = res.Bivector0;
        }
        if (elems::has_bivectorE(elements))
        {
            out.BivectorE = res.BivectorE;
        }
        if (elems::has_trivector(elements))
        {
            out.Trivector = res.Trivector;
        }
        return out;
    }

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
