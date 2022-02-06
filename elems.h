/*
 * This file is part of the Tiny-PGA distribution (https://github.com/sergehog/tiny_pga)
 * Copyright (c) 2022 Sergey Smirnov / Seregium Oy.
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

#ifndef TINY_PGA_ELEMS_H
#define TINY_PGA_ELEMS_H

#include <array>
#include <cstdint>

namespace tiny_pga
{

/// Internal namespace for compile-time operations with multivector elements
namespace elems
{

/// Bitmap of elements in the Multivector
using Elems = std::uint16_t;

/// Element names,  enumerated in some logical order as bit indexes
enum Names : std::uint8_t
{
    // Vector
    E1 = 0U,
    E2 = 1U,
    E3 = 2U,
    E0 = 3U,

    // BivectorE / Quaternion / Rotor
    Scalar = 4U,
    E12 = 5U,
    E31 = 6U,
    E23 = 7U,

    // Bivector0 / Translator
    E01 = 8U,
    E02 = 9U,
    E03 = 10U,
    E0123 = 11U,  // Pseudo Scalar

    // Trivector / Point
    E021 = 12U,
    E013 = 13U,
    E032 = 14U,
    E123 = 15U,

    Amount = 16U
};

/// Bitset for individual mulivector elements
enum class Values : std::uint16_t
{
    // Vector : corresponds to Plane
    kE1 = 1U << uint16_t(E1),
    kE2 = 1U << uint16_t(E2),
    kE3 = 1U << uint16_t(E3),
    kE0 = 1U << uint16_t(E0),

    // Bi-vectorE : corresponds to Quaternion / Rotor
    kScalar = 1U << uint16_t(Scalar),
    kE12 = 1U << uint16_t(E12),
    kE31 = 1U << uint16_t(E31),
    kE23 = 1U << uint16_t(E23),

    // Bi-vector0 : corresponds to Translator
    kE01 = 1U << uint16_t(E01),
    kE02 = 1U << uint16_t(E02),
    kE03 = 1U << uint16_t(E03),
    kE0123 = 1U << uint16_t(E0123),  // Pseudo Scalar

    // Trivector : corresponds to Point
    kE021 = 1U << uint16_t(E021),
    kE013 = 1U << uint16_t(E013),
    kE032 = 1U << uint16_t(E032),
    kE123 = 1U << uint16_t(E123),
};

/// Check if bitset of Elements has particular element in it
template <Names name>
constexpr bool has_elem(const Elems elems)
{
    return (elems & static_cast<Elems>(1U << uint16_t(name))) > 0;
}

/// Index of the particular element, in the multivector array
template <Names name>
constexpr std::size_t index(const Elems elems)
{
    std::size_t index{};
    for (std::uint8_t i = 0U; i < std::uint8_t(name); i++)
    {
        if ((elems & (1U << uint16_t(i))) > 0)
        {
            index++;
        }
    }
    return index;
}

constexpr std::array<std::size_t, elems::Names::Amount> indexes(const elems::Elems elements)
{
    using Indexes = std::array<std::size_t, elems::Names::Amount>;
    Indexes indexes{index<E1>(elements),
                    index<E2>(elements),
                    index<E3>(elements),
                    index<E0>(elements),

                    index<Scalar>(elements),
                    index<E12>(elements),
                    index<E31>(elements),
                    index<E23>(elements),

                    index<E01>(elements),
                    index<E02>(elements),
                    index<E03>(elements),
                    index<E0123>(elements),

                    index<E021>(elements),
                    index<E013>(elements),
                    index<E032>(elements),
                    index<E123>(elements)};
    return indexes;
}

/// Number of set elements (multivector array size)
constexpr std::size_t count(Elems elems)
{
    std::size_t count{};
    for (std::uint8_t i = 0U; i < std::uint8_t(Names::Amount); i++)
    {
        if ((elems & (1U << uint16_t(i))) > 0)
        {
            count++;
        }
    }
    return count;
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
    out_elements |= scalar ? static_cast<Elems>(Values::kScalar) : 0U;
    out_elements |= e0 ? static_cast<Elems>(Values::kE0) : 0U;
    out_elements |= e1 ? static_cast<Elems>(Values::kE1) : 0U;
    out_elements |= e2 ? static_cast<Elems>(Values::kE2) : 0U;
    out_elements |= e3 ? static_cast<Elems>(Values::kE3) : 0U;
    out_elements |= e01 ? static_cast<Elems>(Values::kE01) : 0U;
    out_elements |= e02 ? static_cast<Elems>(Values::kE02) : 0U;
    out_elements |= e03 ? static_cast<Elems>(Values::kE03) : 0U;
    out_elements |= e12 ? static_cast<Elems>(Values::kE12) : 0U;
    out_elements |= e31 ? static_cast<Elems>(Values::kE31) : 0U;
    out_elements |= e23 ? static_cast<Elems>(Values::kE23) : 0U;
    out_elements |= e021 ? static_cast<Elems>(Values::kE021) : 0U;
    out_elements |= e013 ? static_cast<Elems>(Values::kE013) : 0U;
    out_elements |= e032 ? static_cast<Elems>(Values::kE032) : 0U;
    out_elements |= e123 ? static_cast<Elems>(Values::kE123) : 0U;
    out_elements |= e0123 ? static_cast<Elems>(Values::kE0123) : 0U;
    return out_elements;
}

/// Some shorter aliases for frequently used functions
constexpr auto has_scalar = has_elem<Scalar>;
constexpr auto has_e0 = has_elem<E0>;
constexpr auto has_e1 = has_elem<E1>;
constexpr auto has_e2 = has_elem<E2>;
constexpr auto has_e3 = has_elem<E3>;
constexpr auto has_e12 = has_elem<E12>;
constexpr auto has_e23 = has_elem<E23>;
constexpr auto has_e31 = has_elem<E31>;
constexpr auto has_e01 = has_elem<E01>;
constexpr auto has_e02 = has_elem<E02>;
constexpr auto has_e03 = has_elem<E03>;
constexpr auto has_e021 = has_elem<E021>;
constexpr auto has_e013 = has_elem<E013>;
constexpr auto has_e032 = has_elem<E032>;
constexpr auto has_e123 = has_elem<E123>;
constexpr auto has_e0123 = has_elem<E0123>;

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

constexpr Elems inner_product(const Elems elems1, const Elems elems2)
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
                    (has_e2(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e2(elems2)) ||
                    (has_e3(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e3(elems2)) ||
                    (has_e23(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e23(elems2));

    const bool e2 = (has_scalar(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e1(elems2)) ||
                    (has_e3(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e3(elems2)) ||
                    (has_e31(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e31(elems2));

    const bool e3 = (has_scalar(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e1(elems2)) ||
                    (has_e2(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e2(elems2)) ||
                    (has_e12(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e12(elems2));

    const bool e01 = (has_scalar(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_scalar(elems2)) ||
                     (has_e2(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e2(elems2)) ||
                     (has_e3(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e3(elems2)) ||
                     (has_e23(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e23(elems2));

    const bool e02 = (has_scalar(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e1(elems2)) ||
                     (has_e3(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e3(elems2)) ||
                     (has_e31(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e31(elems2));

    const bool e03 = (has_scalar(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e1(elems2)) ||
                     (has_e2(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e2(elems2)) ||
                     (has_e12(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e12(elems2));

    const bool e12 = (has_scalar(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_scalar(elems2)) ||
                     (has_e3(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e3(elems2));

    const bool e23 = (has_scalar(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e1(elems2));

    const bool e31 = (has_scalar(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_scalar(elems2)) ||
                     (has_e2(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e2(elems2));

    const bool e021 = (has_scalar(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_scalar(elems2)) ||
                      (has_e3(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e3(elems2));

    const bool e013 = (has_scalar(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_scalar(elems2)) ||
                      (has_e2(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e2(elems2));

    const bool e032 = (has_scalar(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_scalar(elems2)) ||
                      (has_e1(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e1(elems2));

    const bool e123 = (has_scalar(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_scalar(elems2));

    const bool e0123 = (has_scalar(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_scalar(elems2));

    return elements(scalar, e0, e1, e2, e3, e01, e02, e03, e12, e31, e23, e021, e013, e032, e123, e0123);
}

constexpr Elems outer_product(const Elems elems1, const Elems elems2)
{
    const bool scalar = (has_scalar(elems1) && has_scalar(elems2));

    const bool e0 = (has_scalar(elems1) && has_e0(elems2)) || (has_e0(elems1) && has_scalar(elems2));

    const bool e1 = (has_scalar(elems1) && has_e1(elems2)) || (has_e1(elems1) && has_scalar(elems2));

    const bool e2 = (has_scalar(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_scalar(elems2));

    const bool e3 = (has_scalar(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_scalar(elems2));

    const bool e01 = (has_scalar(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e1(elems2)) || (has_e1(elems1) && has_e0(elems2));

    const bool e02 = (has_scalar(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_e0(elems2));

    const bool e03 = (has_scalar(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_scalar(elems2)) ||
                     (has_e0(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e0(elems2));

    const bool e12 = (has_scalar(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e2(elems2)) || (has_e2(elems1) && has_e1(elems2));

    const bool e31 = (has_scalar(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_scalar(elems2)) ||
                     (has_e1(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e1(elems2));

    const bool e23 = (has_scalar(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_scalar(elems2)) ||
                     (has_e2(elems1) && has_e3(elems2)) || (has_e3(elems1) && has_e2(elems2));

    const bool e021 = (has_scalar(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e0(elems2)) ||
                      (has_e1(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e1(elems2)) ||
                      (has_e2(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e2(elems2));

    const bool e013 = (has_scalar(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e0(elems2)) ||
                      (has_e1(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e1(elems2)) ||
                      (has_e3(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e3(elems2));

    const bool e032 = (has_scalar(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_scalar(elems2)) ||
                      (has_e0(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e0(elems2)) ||
                      (has_e2(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e2(elems2)) ||
                      (has_e3(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e3(elems2));

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

constexpr Elems regressive_product(const Elems elems1, const Elems elems2)
{
    return 0;
}

constexpr Elems commutator_product(const Elems elems1, const Elems elems2)
{
    return 0;
}

constexpr Elems sandwich_product(const Elems elems1, const Elems elems2)
{
    return 0;
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

constexpr Elems dual(Elems elems)
{
    Elems result{0U};
    result |= has_e0123(elems) ? static_cast<Elems>(elems::Values::kScalar) : 0UL;
    result |= has_e123(elems) ? static_cast<Elems>(elems::Values::kE0) : 0UL;
    result |= has_e032(elems) ? static_cast<Elems>(elems::Values::kE1) : 0UL;
    result |= has_e013(elems) ? static_cast<Elems>(elems::Values::kE2) : 0UL;
    result |= has_e021(elems) ? static_cast<Elems>(elems::Values::kE3) : 0UL;

    result |= has_e23(elems) ? static_cast<Elems>(elems::Values::kE01) : 0UL;
    result |= has_e31(elems) ? static_cast<Elems>(elems::Values::kE02) : 0UL;
    result |= has_e12(elems) ? static_cast<Elems>(elems::Values::kE03) : 0UL;

    result |= has_e01(elems) ? static_cast<Elems>(elems::Values::kE23) : 0UL;
    result |= has_e02(elems) ? static_cast<Elems>(elems::Values::kE31) : 0UL;
    result |= has_e03(elems) ? static_cast<Elems>(elems::Values::kE12) : 0UL;

    result |= has_e0(elems) ? static_cast<Elems>(elems::Values::kE123) : 0UL;
    result |= has_e1(elems) ? static_cast<Elems>(elems::Values::kE032) : 0UL;
    result |= has_e2(elems) ? static_cast<Elems>(elems::Values::kE013) : 0UL;
    result |= has_e3(elems) ? static_cast<Elems>(elems::Values::kE021) : 0UL;

    result |= has_scalar(elems) ? static_cast<Elems>(elems::Values::kE0123) : 0UL;
    return result;
}

constexpr Elems PlaneElems = static_cast<Elems>(elems::Values::kE0) + static_cast<Elems>(elems::Values::kE1) +
                             static_cast<Elems>(elems::Values::kE2) + static_cast<Elems>(elems::Values::kE3);

constexpr Elems ComplexElems = static_cast<Elems>(elems::Values::kScalar) + static_cast<Elems>(elems::Values::kE12);

constexpr Elems LineElems = static_cast<Elems>(elems::Values::kE01) + static_cast<Elems>(elems::Values::kE02) +
                            static_cast<Elems>(elems::Values::kE03) + static_cast<Elems>(elems::Values::kE23) +
                            static_cast<Elems>(elems::Values::kE31) + static_cast<Elems>(elems::Values::kE12);

constexpr Elems PointElems = static_cast<Elems>(elems::Values::kE123) + static_cast<Elems>(elems::Values::kE032) +
                             static_cast<Elems>(elems::Values::kE013) + static_cast<Elems>(elems::Values::kE021);

constexpr Elems RotorElems = static_cast<Elems>(elems::Values::kScalar) + static_cast<Elems>(elems::Values::kE23) +
                             static_cast<Elems>(elems::Values::kE31) + static_cast<Elems>(elems::Values::kE12);

constexpr Elems TranslatorElems = static_cast<Elems>(elems::Values::kE01) + static_cast<Elems>(elems::Values::kE02) +
                                  static_cast<Elems>(elems::Values::kE03) + static_cast<Elems>(elems::Values::kScalar);

constexpr Elems MotorElems = static_cast<Elems>(elems::Values::kScalar) + static_cast<Elems>(elems::Values::kE23) +
                             static_cast<Elems>(elems::Values::kE31) + static_cast<Elems>(elems::Values::kE12) +
                             static_cast<Elems>(elems::Values::kE01) + static_cast<Elems>(elems::Values::kE02) +
                             static_cast<Elems>(elems::Values::kE03) + static_cast<Elems>(elems::Values::kE0123);

}  // namespace elems
}  // namespace tiny_pga

#endif  // TINY_PGA_ELEMS_H
