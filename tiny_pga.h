/*
 * This file is part of the Tiny-PGA distribution (https://github.com/sergehog/tiny_pga)
 * Copyright (c) 2020-2022 Sergey Smirnov / Seregium Oy.
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
// C++ implementation of 3D Projective Geometric Algebra (aka Plane-based Geometric Algebra)
// C++ Templating helps with number of issues:
// * reducing computational complexity of PGA operators via compile-time optimizations
// * reducing memory footprint
// * compile-time checks for assignments correctness / blade matching (works to some extent)
//

#ifndef TINY_PGA_H
#define TINY_PGA_H

#include "elems.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <type_traits>

namespace tiny_pga
{

/// Compile-time optimized (using C++ templating) implementation of 3D PGA Multivector
template <elems::Elems elements, typename ScalarType = float>
struct Multivector;

/// Geometric Product (generic implementation)
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::geometric_product(first_elements, second_elements), ScalarType> geometric_product(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Inner (dot) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::inner_product(first_elements, second_elements), ScalarType> inner_product(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Outer (meet) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::outer_product(first_elements, second_elements), ScalarType> outer_product(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Regressive (join) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::regressive_product(first_elements, second_elements), ScalarType> regressive_product(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Commutator (cross) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::commutator_product(first_elements, second_elements), ScalarType> commutator_product(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Sandwich (transform) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<first_elements, ScalarType> sandwich_product(const Multivector<first_elements, ScalarType>&,
                                                         const Multivector<second_elements, ScalarType>&);

/// Addition
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> addition(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Subtraction
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> subtraction(
    const Multivector<first_elements, ScalarType>&,
    const Multivector<second_elements, ScalarType>&);

/// Reverse
template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> reverse(const Multivector<elements, ScalarType>&);

/// Dual
template <elems::Elems elements, typename ScalarType>
Multivector<elems::dual(elements), ScalarType> dual(const Multivector<elements, ScalarType>&);

/// Clifford Conjugation
template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> conjugate(const Multivector<elements, ScalarType>& object);

template <elems::Elems elements, typename ScalarType>
struct Multivector
{
    using Type = ScalarType;
    /// Exposing elements outside, as static const value
    static constexpr elems::Elems Elements = elements;

    /// values of Multivector elements
    std::array<ScalarType, elems::count(elements)> values{};

    Multivector(std::initializer_list<ScalarType> list) { std::copy(list.begin(), list.end(), values.begin()); }

    /// explicit cast from other Multivector type
    template <elems::Elems other_elements, typename OtherScalarType>
    explicit Multivector(const Multivector<other_elements, OtherScalarType>& other)
    {
        using namespace elems;
        if constexpr (has_elem<Scalar>(elements) && has_elem<Scalar>(other_elements))
        {
            value<Scalar>() = static_cast<ScalarType>(other.template value<Scalar>());
        }
        if constexpr (has_elem<E0>(elements) && has_elem<E0>(other_elements))
        {
            value<E0>() = static_cast<ScalarType>(other.template value<E0>());
        }
        if constexpr (has_elem<E1>(elements) && has_elem<E1>(other_elements))
        {
            value<E1>() = static_cast<ScalarType>(other.template value<E1>());
        }
        if constexpr (has_elem<E2>(elements) && has_elem<E2>(other_elements))
        {
            value<E2>() = static_cast<ScalarType>(other.template value<E2>());
        }
        if constexpr (has_elem<E3>(elements) && has_elem<E3>(other_elements))
        {
            value<E3>() = static_cast<ScalarType>(other.template value<E3>());
        }
        if constexpr (has_elem<E01>(elements) && has_elem<E01>(other_elements))
        {
            value<E01>() = static_cast<ScalarType>(other.template value<E01>());
        }
        if constexpr (has_elem<E02>(elements) && has_elem<E02>(other_elements))
        {
            value<E02>() = static_cast<ScalarType>(other.template value<E02>());
        }
        if constexpr (has_elem<E03>(elements) && has_elem<E03>(other_elements))
        {
            value<E03>() = static_cast<ScalarType>(other.template value<E03>());
        }
        if constexpr (has_elem<E12>(elements) && has_elem<E12>(other_elements))
        {
            value<E12>() = static_cast<ScalarType>(other.template value<E12>());
        }
        if constexpr (has_elem<E23>(elements) && has_elem<E23>(other_elements))
        {
            value<E23>() = static_cast<ScalarType>(other.template value<E23>());
        }
        if constexpr (has_elem<E31>(elements) && has_elem<E31>(other_elements))
        {
            value<E31>() = static_cast<ScalarType>(other.template value<E31>());
        }
        if constexpr (has_elem<E032>(elements) && has_elem<E032>(other_elements))
        {
            value<E032>() = static_cast<ScalarType>(other.template value<E032>());
        }
        if constexpr (has_elem<E013>(elements) && has_elem<E013>(other_elements))
        {
            value<E013>() = static_cast<ScalarType>(other.template value<E013>());
        }
        if constexpr (has_elem<E021>(elements) && has_elem<E021>(other_elements))
        {
            value<E021>() = static_cast<ScalarType>(other.template value<E021>());
        }
        if constexpr (has_elem<E123>(elements) && has_elem<E123>(other_elements))
        {
            value<E123>() = static_cast<ScalarType>(other.template value<E123>());
        }
        if constexpr (has_elem<E0123>(elements) && has_elem<E0123>(other_elements))
        {
            value<E0123>() = static_cast<ScalarType>(other.template value<E0123>());
        }
    }

    /// indexes
    // static constexpr std::array<std::size_t, elems::Names::Amount> indexes{elems::indexes(elements)};

    /// compile-time getter/setter for particular element
    template <elems::Names elem>
    ScalarType& value()
    {
        static_assert(elems::has_elem<elem>(elements), "Multivector has no such element");
        return values[elems::index<elem>(elements)];
    }

    /// compile-time getter for particular element
    template <elems::Names elem>
    ScalarType value() const
    {
        static_assert(elems::has_elem<elem>(elements), "Multivector has no such element");
        return values.at(elems::index<elem>(elements));
    }

    /// runtime access (getter/setter) to elements
    ScalarType& operator[](elems::Names elem)
    {
        size_t index{};
        for (std::uint8_t i = 0U; i < std::uint8_t(elem); i++)
        {
            if ((elements & (1U << uint16_t(i))) > 0)
            {
                index++;
            }
        }
        return values[index];
    }

    /// runtime access (getter) to elements
    ScalarType operator[](elems::Names elem) const
    {
        size_t index{};
        for (std::uint8_t i = 0U; i < std::uint8_t(elem); i++)
        {
            if ((elements & (1U << uint16_t(i))) > 0)
            {
                index++;
            }
        }
        return values[index];
    }

    /// Geometric product operator
    template <elems::Elems other_elements>
    Multivector<elems::geometric_product(elements, other_elements), ScalarType> operator*(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::geometric_product(*this, other);
    }

    /// Inner product operator (also known as Dot-product)
    template <elems::Elems other_elements>
    Multivector<elems::inner_product(elements, other_elements), ScalarType> operator|(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::inner_product(*this, other);
    }

    /// Alias for Inner product (Dot-product)
    template <elems::Elems other_elements>
    Multivector<elems::inner_product(elements, other_elements), ScalarType> dot(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::inner_product(*this, other);
    }

    /// Outer product operator (also known as Wedge- or Meet- product)
    template <elems::Elems other_elements>
    Multivector<elems::outer_product(elements, other_elements), ScalarType> operator^(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::outer_product(*this, other);
    }

    /// Alias for Outer product operator
    template <elems::Elems other_elements>
    Multivector<elems::outer_product(elements, other_elements), ScalarType> meet(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::outer_product(*this, other);
    }

    /// Regressive product operator (also known as Vee- or Join-product)
    template <elems::Elems other_elements>
    Multivector<elems::regressive_product(elements, other_elements), ScalarType> operator&(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::regressive_product(*this, other);
    }

    /// Alias for Regressive product operator
    template <elems::Elems other_elements>
    Multivector<elems::regressive_product(elements, other_elements), ScalarType> join(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::regressive_product(*this, other);
    }

    /// Multivector addition
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> operator+(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::addition(*this, other);
    }

    /// Alias for addition
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> add(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::addition(*this, other);
    }

    /// Multivector subtraction
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> operator-(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::subtraction(*this, other);
    }

    /// Unary minus operator
    Multivector<elements, ScalarType> operator-() const
    {
        Multivector<0, ScalarType> nothing{};
        return tiny_pga::subtraction(nothing, *this);
    }

    /// Alias for subtraction
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> sub(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::subtraction(*this, other);
    }

    /// Reverse operator
    Multivector<elements, ScalarType> operator~() const { return tiny_pga::reverse(*this); }

    /// Reverse function
    Multivector<elements, ScalarType> reverse() const { return tiny_pga::reverse(*this); }

    /// Duality operator
    Multivector<elems::dual(elements), ScalarType> operator!() const { return tiny_pga::dual(*this); }

    /// Duality function
    Multivector<elems::dual(elements), ScalarType> dual() const { return tiny_pga::dual(*this); }

    /// sandwich product with Rotor
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::RotorElems, ScalarType>& other) const
    {
        return tiny_pga::sandwich_product(*this, other);
    }

    /// sandwich product with Translator
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::TranslatorElems, ScalarType>& other) const
    {
        return tiny_pga::sandwich_product(*this, other);
    }

    /// sandwich product with Motor
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::MotorElems, ScalarType>& other) const
    {
        return tiny_pga::sandwich_product(*this, other);
    }

    Multivector<elements, ScalarType> Normalized() const { return *this * ScalarType(1.F / norm(*this)); }
};

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::geometric_product(first_elements, second_elements), ScalarType> operator*(
    const Multivector<first_elements, ScalarType>& a,
    const Multivector<second_elements, ScalarType>& b)
{
    return tiny_pga::geometric_product(a, b);
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, second_elements), ScalarType> operator+(
    const Multivector<first_elements, ScalarType>& a,
    const Multivector<second_elements, ScalarType>& b)
{
    return tiny_pga::addition(a, b);
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, second_elements), ScalarType> operator-(
    const Multivector<first_elements, ScalarType>& a,
    const Multivector<second_elements, ScalarType>& b)
{
    return tiny_pga::subtraction(a, b);
}

#ifdef TRY_AVOIDING_MACCOS
// instead of ELEM_MULTIPLY / ELEM_BOTH_MULTIPLY maxcros we may try something like that
// namespace internal
//{
// template <elems::Elems first_elements,
//          elems::Elems second_elements,
//          elems::Names elem1,
//          elems::Names elem2,
//          typename ScalarType>
// ScalarType elemMultiply(const Multivector<first_elements, ScalarType>& first,
//                        const Multivector<second_elements, ScalarType>& second)
//{
//    if constexpr (elems::has_elem<elem1>(first_elements) && elems::has_elem<elem2>(second_elements))
//    {
//        return  first.template value<elem1>() * second.template value<elem2>();
//    }
//    return 0;
//}
// out.template value<Scalar>() += internal::elemMultiply<first_elements, second_elements, Scalar, Scalar,
// ScalarType>(first, second); out.template value<Scalar>() += internal::elemMultiply<first_elements, second_elements,
// E1, E1, ScalarType>(first, second); out.template value<Scalar>() += internal::elemMultiply<first_elements,
// second_elements, E2, E2, ScalarType>(first, second); out.template value<Scalar>() +=
// internal::elemMultiply<first_elements, second_elements, E3, E3, ScalarType>(first, second); out.template
// value<Scalar>() -= internal::elemMultiply<first_elements, second_elements, E12, E12, ScalarType>(first, second);
// out.template value<Scalar>() -= internal::elemMultiply<first_elements, second_elements, E31, E31, ScalarType>(first,
// second); out.template value<Scalar>() -= internal::elemMultiply<first_elements, second_elements, E23, E23,
// ScalarType>(first, second); out.template value<Scalar>() -= internal::elemMultiply<first_elements, second_elements,
// E123, E123, ScalarType>(first, second); } //namespace internal
#endif

/// Geometric Product operation
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::geometric_product(first_elements, second_elements), ScalarType> geometric_product(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    using namespace elems;
    constexpr Elems out_elems = geometric_product(first_elements, second_elements);
    Multivector<out_elems, ScalarType> out{};

#define ELEM_MULTIPLY(elem_out, elem1, elem2, sign1)                                                         \
    if constexpr (has_elem<elem1>(first_elements) && has_elem<elem2>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign1 first.template value<elem1>() * second.template value<elem2>(); \
    }

#define ELEM_BOTH_MULTIPLY(elem_out, elem1, elem2, sign1, sign2)                                             \
    if constexpr (has_elem<elem1>(first_elements) && has_elem<elem2>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign1 first.template value<elem1>() * second.template value<elem2>(); \
    }                                                                                                        \
    if constexpr (has_elem<elem2>(first_elements) && has_elem<elem1>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign2 first.template value<elem2>() * second.template value<elem1>(); \
    }

    if constexpr (has_scalar(out_elems))
    {
        ELEM_MULTIPLY(Scalar, Scalar, Scalar, +=);
        ELEM_MULTIPLY(Scalar, E1, E1, +=);
        ELEM_MULTIPLY(Scalar, E2, E2, +=);
        ELEM_MULTIPLY(Scalar, E3, E3, +=);
        ELEM_MULTIPLY(Scalar, E12, E12, -=);
        ELEM_MULTIPLY(Scalar, E31, E31, -=);
        ELEM_MULTIPLY(Scalar, E23, E23, -=);
        ELEM_MULTIPLY(Scalar, E123, E123, -=);
    }
#undef ELEM_MULTIPLY

    if constexpr (has_e0(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0, Scalar, E0, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E1, E01, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E2, E02, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E3, E03, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E12, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E31, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E23, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E123, E0123, +=, -=);
    }

    if constexpr (has_e1(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E1, Scalar, E1, +=, +=);
        ELEM_BOTH_MULTIPLY(E1, E2, E12, -=, +=);
        ELEM_BOTH_MULTIPLY(E1, E3, E31, +=, -=);
        ELEM_BOTH_MULTIPLY(E1, E23, E123, -=, -=);
    }

    if constexpr (has_e2(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E2, Scalar, E2, +=, +=);
        ELEM_BOTH_MULTIPLY(E2, E1, E12, +=, -=);
        ELEM_BOTH_MULTIPLY(E2, E3, E23, -=, +=);
        ELEM_BOTH_MULTIPLY(E2, E31, E123, -=, -=);
    }

    if constexpr (has_e3(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E3, Scalar, E3, +=, +=);
        ELEM_BOTH_MULTIPLY(E3, E1, E31, -=, +=);
        ELEM_BOTH_MULTIPLY(E3, E2, E23, +=, -=);
        ELEM_BOTH_MULTIPLY(E3, E12, E123, -=, -=);
    }

    if constexpr (has_e01(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E01, Scalar, E01, +=, +=);
        ELEM_BOTH_MULTIPLY(E01, E0, E1, +=, -=);
        ELEM_BOTH_MULTIPLY(E01, E2, E021, -=, -=);
        ELEM_BOTH_MULTIPLY(E01, E3, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E01, E02, E12, -=, +=);
        ELEM_BOTH_MULTIPLY(E01, E03, E31, +=, -=);
        ELEM_BOTH_MULTIPLY(E01, E23, E0123, -=, -=);
        ELEM_BOTH_MULTIPLY(E01, E032, E123, +=, -=);
    }

    if constexpr (has_e02(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E02, Scalar, E02, +=, +=);
        ELEM_BOTH_MULTIPLY(E02, E0, E2, +=, -=);
        ELEM_BOTH_MULTIPLY(E02, E1, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E02, E3, E032, -=, -=);
        ELEM_BOTH_MULTIPLY(E02, E01, E12, +=, -=);
        ELEM_BOTH_MULTIPLY(E02, E03, E23, -=, +=);
        ELEM_BOTH_MULTIPLY(E02, E31, E0123, -=, -=);
        ELEM_BOTH_MULTIPLY(E02, E013, E123, +=, -=);
    }

    if constexpr (has_e03(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E03, Scalar, E03, +=, +=);
        ELEM_BOTH_MULTIPLY(E03, E0, E3, +=, -=);
        ELEM_BOTH_MULTIPLY(E03, E1, E013, -=, -=);
        ELEM_BOTH_MULTIPLY(E03, E2, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E03, E01, E31, -=, +=);
        ELEM_BOTH_MULTIPLY(E03, E02, E23, +=, -=);
        ELEM_BOTH_MULTIPLY(E03, E12, E0123, -=, -=);
        ELEM_BOTH_MULTIPLY(E03, E021, E123, +=, -=);
    }

    if constexpr (has_e12(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E12, Scalar, E12, +=, +=);
        ELEM_BOTH_MULTIPLY(E12, E1, E2, +=, -=);
        ELEM_BOTH_MULTIPLY(E12, E3, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E12, E31, E23, +=, -=);
    }

    if constexpr (has_e31(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E31, Scalar, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E31, E3, E1, +=, -=);
        ELEM_BOTH_MULTIPLY(E31, E2, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E31, E12, E23, -=, +=);
    }

    if constexpr (has_e23(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E23, Scalar, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E23, E2, E3, +=, -=);
        ELEM_BOTH_MULTIPLY(E23, E1, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E23, E12, E31, +=, -=);
    }

    if constexpr (has_e021(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E021, Scalar, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E021, E0, E12, +=, -=);
        ELEM_BOTH_MULTIPLY(E021, E1, E02, +=, +=);
        ELEM_BOTH_MULTIPLY(E021, E2, E01, -=, -=);
        ELEM_BOTH_MULTIPLY(E021, E3, E0123, +=, -=);
        ELEM_BOTH_MULTIPLY(E021, E03, E123, +=, -=);
        ELEM_BOTH_MULTIPLY(E021, E31, E032, +=, -=);
        ELEM_BOTH_MULTIPLY(E021, E23, E013, -=, +=);
    }

    if constexpr (has_e013(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E013, Scalar, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E013, E0, E31, +=, -=);
        ELEM_BOTH_MULTIPLY(E013, E1, E03, +=, -=);
        ELEM_BOTH_MULTIPLY(E013, E3, E01, +=, +=);
        ELEM_BOTH_MULTIPLY(E013, E2, E0123, +=, -=);
        ELEM_BOTH_MULTIPLY(E013, E02, E123, +=, -=);
        ELEM_BOTH_MULTIPLY(E013, E12, E032, -=, +=);
        ELEM_BOTH_MULTIPLY(E013, E23, E021, -=, +=);
    }

    if constexpr (has_e032(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E032, Scalar, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E032, E0, E23, -=, -=);
        ELEM_BOTH_MULTIPLY(E032, E2, E03, +=, +=);
        ELEM_BOTH_MULTIPLY(E032, E3, E02, -=, -=);
        ELEM_BOTH_MULTIPLY(E032, E1, E0123, +=, -=);
        ELEM_BOTH_MULTIPLY(E032, E01, E123, -=, +=);
        ELEM_BOTH_MULTIPLY(E032, E12, E013, +=, -=);
        ELEM_BOTH_MULTIPLY(E032, E31, E021, -=, +=);
    }

    if constexpr (has_e123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E123, Scalar, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E1, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E2, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E3, E12, +=, +=);
    }

    if constexpr (has_e0123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0123, Scalar, E0123, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E0, E123, +=, -=);
        ELEM_BOTH_MULTIPLY(E0123, E1, E032, +=, -=);
        ELEM_BOTH_MULTIPLY(E0123, E2, E013, +=, -=);
        ELEM_BOTH_MULTIPLY(E0123, E3, E021, +=, -=);
        ELEM_BOTH_MULTIPLY(E0123, E01, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E02, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E03, E12, +=, +=);
    }
#undef ELEM_BOTH_MULTIPLY
    return out;
}

/// Inner Product operation
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::inner_product(first_elements, second_elements), ScalarType> inner_product(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    using namespace elems;
    constexpr Elems out_elems = inner_product(first_elements, second_elements);
    Multivector<out_elems, ScalarType> out{};

#define ELEM_MULTIPLY(elem_out, elem1, elem2, sign1)                                                         \
    if constexpr (has_elem<elem1>(first_elements) && has_elem<elem2>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign1 first.template value<elem1>() * second.template value<elem2>(); \
    }

#define ELEM_BOTH_MULTIPLY(elem_out, elem1, elem2, sign1, sign2)                                             \
    if constexpr (has_elem<elem1>(first_elements) && has_elem<elem2>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign1 first.template value<elem1>() * second.template value<elem2>(); \
    }                                                                                                        \
    if constexpr (has_elem<elem2>(first_elements) && has_elem<elem1>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign2 first.template value<elem2>() * second.template value<elem1>(); \
    }

    if constexpr (has_scalar(out_elems))
    {
        ELEM_MULTIPLY(Scalar, Scalar, Scalar, +=);
        ELEM_MULTIPLY(Scalar, E1, E1, +=);
        ELEM_MULTIPLY(Scalar, E2, E2, +=);
        ELEM_MULTIPLY(Scalar, E3, E3, +=);
        ELEM_MULTIPLY(Scalar, E12, E12, -=);
        ELEM_MULTIPLY(Scalar, E31, E31, -=);
        ELEM_MULTIPLY(Scalar, E23, E23, -=);
        ELEM_MULTIPLY(Scalar, E123, E123, -=);
    }
#undef ELEM_MULTIPLY

    if constexpr (has_e0(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0, Scalar, E0, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E1, E01, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E2, E02, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E3, E03, -=, +=);
        ELEM_BOTH_MULTIPLY(E0, E12, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E31, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E23, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E0, E123, E0123, +=, -=);
    }

    if constexpr (has_e1(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E1, Scalar, E1, +=, +=);
        ELEM_BOTH_MULTIPLY(E1, E2, E12, -=, +=);
        ELEM_BOTH_MULTIPLY(E1, E3, E31, +=, -=);
        ELEM_BOTH_MULTIPLY(E1, E23, E123, -=, -=);
    }

    if constexpr (has_e2(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E2, Scalar, E2, +=, +=);
        ELEM_BOTH_MULTIPLY(E2, E1, E12, +=, -=);
        ELEM_BOTH_MULTIPLY(E2, E3, E23, -=, +=);
        ELEM_BOTH_MULTIPLY(E2, E31, E123, -=, -=);
    }

    if constexpr (has_e3(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E3, Scalar, E3, +=, +=);
        ELEM_BOTH_MULTIPLY(E3, E1, E31, -=, +=);
        ELEM_BOTH_MULTIPLY(E3, E2, E23, +=, -=);
        ELEM_BOTH_MULTIPLY(E3, E12, E123, -=, -=);
    }

    if constexpr (has_e01(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E01, Scalar, E01, +=, +=);
        ELEM_BOTH_MULTIPLY(E01, E2, E021, -=, -=);
        ELEM_BOTH_MULTIPLY(E01, E3, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E01, E032, E123, +=, -=);
    }

    if constexpr (has_e02(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E02, Scalar, E02, +=, +=);
        ELEM_BOTH_MULTIPLY(E02, E1, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E02, E3, E032, -=, -=);
        ELEM_BOTH_MULTIPLY(E02, E31, E0123, -=, -=);
    }

    if constexpr (has_e03(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E03, Scalar, E03, +=, +=);
        ELEM_BOTH_MULTIPLY(E03, E1, E013, -=, -=);
        ELEM_BOTH_MULTIPLY(E03, E2, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E03, E021, E123, +=, -=);
    }

    if constexpr (has_e12(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E12, Scalar, E12, +=, +=);
        ELEM_BOTH_MULTIPLY(E12, E3, E123, +=, +=);
    }

    if constexpr (has_e31(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E31, Scalar, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E31, E2, E123, +=, +=);
    }

    if constexpr (has_e23(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E23, Scalar, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E23, E1, E123, +=, +=);
    }

    if constexpr (has_e021(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E021, Scalar, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E021, E3, E0123, +=, -=);
    }

    if constexpr (has_e013(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E013, Scalar, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E013, E2, E0123, +=, -=);
    }

    if constexpr (has_e032(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E032, Scalar, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E032, E1, E0123, +=, -=);
    }

    if constexpr (has_e123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E123, Scalar, E123, +=, +=);
    }

    if constexpr (has_e0123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0123, Scalar, E0123, +=, +=);
    }
#undef ELEM_BOTH_MULTIPLY
    return out;
}

/// Outer Product operation
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::outer_product(first_elements, second_elements), ScalarType> outer_product(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    using namespace elems;
    constexpr Elems out_elems = outer_product(first_elements, second_elements);
    Multivector<out_elems, ScalarType> out{};

#define ELEM_BOTH_MULTIPLY(elem_out, elem1, elem2, sign1, sign2)                                             \
    if constexpr (has_elem<elem1>(first_elements) && has_elem<elem2>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign1 first.template value<elem1>() * second.template value<elem2>(); \
    }                                                                                                        \
    if constexpr (has_elem<elem2>(first_elements) && has_elem<elem1>(second_elements))                       \
    {                                                                                                        \
        out.template value<elem_out>() sign2 first.template value<elem2>() * second.template value<elem1>(); \
    }

    if constexpr (has_scalar(out_elems))
    {
        out.template value<Scalar>() += first.template value<Scalar>() + second.template value<Scalar>();
    }

    if constexpr (has_e0(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0, Scalar, E0, +=, +=);
    }

    if constexpr (has_e1(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E1, Scalar, E1, +=, +=);
    }

    if constexpr (has_e2(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E2, Scalar, E2, +=, +=);
    }

    if constexpr (has_e3(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E3, Scalar, E3, +=, +=);
    }

    if constexpr (has_e01(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E01, Scalar, E01, +=, +=);
        ELEM_BOTH_MULTIPLY(E01, E0, E1, +=, -=);
    }

    if constexpr (has_e02(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E02, Scalar, E02, +=, +=);
        ELEM_BOTH_MULTIPLY(E02, E0, E2, +=, -=);
    }

    if constexpr (has_e03(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E03, Scalar, E03, +=, +=);
        ELEM_BOTH_MULTIPLY(E03, E0, E3, +=, -=);
    }

    if constexpr (has_e12(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E12, Scalar, E12, +=, +=);
        ELEM_BOTH_MULTIPLY(E12, E1, E2, +=, -=);
    }

    if constexpr (has_e31(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E31, Scalar, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E31, E3, E1, +=, -=);
    }

    if constexpr (has_e23(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E23, Scalar, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E23, E2, E3, +=, -=);
    }

    if constexpr (has_e021(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E021, Scalar, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E021, E0, E12, -=, +=);
        ELEM_BOTH_MULTIPLY(E021, E1, E02, +=, -=);
        ELEM_BOTH_MULTIPLY(E021, E2, E01, -=, +=);
    }

    if constexpr (has_e013(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E013, Scalar, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E013, E0, E31, -=, +=);
        ELEM_BOTH_MULTIPLY(E013, E1, E03, -=, +=);
        ELEM_BOTH_MULTIPLY(E013, E3, E01, +=, -=);
    }

    if constexpr (has_e032(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E032, Scalar, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E032, E0, E23, -=, +=);
        ELEM_BOTH_MULTIPLY(E032, E2, E03, +=, -=);
        ELEM_BOTH_MULTIPLY(E032, E3, E02, -=, +=);
    }

    if constexpr (has_e123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E123, Scalar, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E1, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E2, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E123, E3, E12, +=, +=);
    }

    if constexpr (has_e0123(out_elems))
    {
        ELEM_BOTH_MULTIPLY(E0123, Scalar, E0123, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E0, E123, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E1, E032, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E2, E013, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E3, E021, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E01, E23, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E02, E31, +=, +=);
        ELEM_BOTH_MULTIPLY(E0123, E03, E12, +=, +=);
    }
#undef ELEM_BOTH_MULTIPLY
    return out;
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> addition(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    constexpr elems::Elems out_elems = elems::addition(first_elements, second_elements);
    Multivector<out_elems, ScalarType> out{};

#define ELEM_ADD(name)                                                   \
    if constexpr (elems::has_elem<name>(out_elems))                      \
    {                                                                    \
        if constexpr (elems::has_elem<name>(first_elements))             \
        {                                                                \
            out.template value<name>() = first.template value<name>();   \
        }                                                                \
        if constexpr (elems::has_elem<name>(second_elements))            \
        {                                                                \
            out.template value<name>() += second.template value<name>(); \
        }                                                                \
    }

    ELEM_ADD(elems::Scalar);
    ELEM_ADD(elems::E0);
    ELEM_ADD(elems::E1);
    ELEM_ADD(elems::E2);
    ELEM_ADD(elems::E3);
    ELEM_ADD(elems::E01);
    ELEM_ADD(elems::E02);
    ELEM_ADD(elems::E03);
    ELEM_ADD(elems::E12);
    ELEM_ADD(elems::E31);
    ELEM_ADD(elems::E23);
    ELEM_ADD(elems::E021);
    ELEM_ADD(elems::E013);
    ELEM_ADD(elems::E032);
    ELEM_ADD(elems::E123);
    ELEM_ADD(elems::E0123);
    return out;
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> subtraction(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    constexpr elems::Elems out_elems = elems::addition(first_elements, second_elements);
    Multivector<out_elems, ScalarType> out{};

#define ELEM_SUB(name)                                                   \
    if constexpr (elems::has_elem<name>(out_elems))                      \
    {                                                                    \
        if constexpr (elems::has_elem<name>(first_elements))             \
        {                                                                \
            out.template value<name>() = first.template value<name>();   \
        }                                                                \
        if constexpr (elems::has_elem<name>(second_elements))            \
        {                                                                \
            out.template value<name>() -= second.template value<name>(); \
        }                                                                \
    }

    ELEM_SUB(elems::Scalar);
    ELEM_SUB(elems::E0);
    ELEM_SUB(elems::E1);
    ELEM_SUB(elems::E2);
    ELEM_SUB(elems::E3);
    ELEM_SUB(elems::E01);
    ELEM_SUB(elems::E02);
    ELEM_SUB(elems::E03);
    ELEM_SUB(elems::E12);
    ELEM_SUB(elems::E31);
    ELEM_SUB(elems::E23);
    ELEM_SUB(elems::E021);
    ELEM_SUB(elems::E013);
    ELEM_SUB(elems::E032);
    ELEM_SUB(elems::E123);
    ELEM_SUB(elems::E0123);
    return out;
}

template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> reverse(const Multivector<elements, ScalarType>& other)
{
    using namespace elems;
    Multivector<elements, ScalarType> out{};
    if constexpr (has_scalar(elements))
    {
        out.template value<Scalar>() = other.template value<Scalar>();
    }
    if constexpr (has_e0(elements))
    {
        out.template value<E0>() = other.template value<E0>();
    }
    if constexpr (has_e1(elements))
    {
        out.template value<E1>() = other.template value<E1>();
    }
    if constexpr (has_e2(elements))
    {
        out.template value<E2>() = other.template value<E2>();
    }
    if constexpr (has_e3(elements))
    {
        out.template value<E3>() = other.template value<E3>();
    }
    if constexpr (has_e01(elements))
    {
        out.template value<E01>() = -other.template value<E01>();
    }
    if constexpr (has_e02(elements))
    {
        out.template value<E02>() = -other.template value<E02>();
    }
    if constexpr (has_e03(elements))
    {
        out.template value<E03>() = -other.template value<E03>();
    }
    if constexpr (has_e12(elements))
    {
        out.template value<E12>() = -other.template value<E12>();
    }
    if constexpr (has_e31(elements))
    {
        out.template value<E31>() = -other.template value<E31>();
    }
    if constexpr (has_e23(elements))
    {
        out.template value<E23>() = -other.template value<E23>();
    }
    if constexpr (has_e021(elements))
    {
        out.template value<E021>() = -other.template value<E021>();
    }
    if constexpr (has_e013(elements))
    {
        out.template value<E013>() = -other.template value<E013>();
    }
    if constexpr (has_e032(elements))
    {
        out.template value<E032>() = -other.template value<E032>();
    }
    if constexpr (has_e123(elements))
    {
        out.template value<E123>() = -other.template value<E123>();
    }
    if constexpr (has_e0123(elements))
    {
        out.template value<E0123>() = other.template value<E0123>();
    }
    return out;
}

template <elems::Elems elements, typename ScalarType>
Multivector<elems::dual(elements), ScalarType> dual(const Multivector<elements, ScalarType>& other)
{
    using namespace elems;
    Multivector<elems::dual(elements), ScalarType> out{};
    if constexpr (elems::has_e0123(elements))
    {
        out.template value<Scalar>() = other.template value<E0123>();
    }
    if constexpr (elems::has_e123(elements))
    {
        out.template value<E0>() = other.template value<E123>();
    }
    if constexpr (elems::has_e032(elements))
    {
        out.template value<E1>() = other.template value<E032>();
    }
    if constexpr (elems::has_e013(elements))
    {
        out.template value<E2>() = other.template value<E013>();
    }
    if constexpr (elems::has_e021(elements))
    {
        out.template value<E3>() = other.template value<E021>();
    }
    if constexpr (elems::has_e23(elements))
    {
        out.template value<E01>() = other.template value<E23>();
    }
    if constexpr (elems::has_e31(elements))
    {
        out.template value<E02>() = other.template value<E31>();
    }
    if constexpr (elems::has_e12(elements))
    {
        out.template value<E03>() = other.template value<E12>();
    }
    if constexpr (elems::has_e01(elements))
    {
        out.template value<E23>() = other.template value<E01>();
    }
    if constexpr (elems::has_e02(elements))
    {
        out.template value<E31>() = other.template value<E02>();
    }
    if constexpr (elems::has_e03(elements))
    {
        out.template value<E12>() = other.template value<E03>();
    }
    if constexpr (elems::has_e0(elements))
    {
        out.template value<E123>() = other.template value<E0>();
    }
    if constexpr (elems::has_e1(elements))
    {
        out.template value<E032>() = other.template value<E1>();
    }
    if constexpr (elems::has_e2(elements))
    {
        out.template value<E013>() = other.template value<E2>();
    }
    if constexpr (elems::has_e3(elements))
    {
        out.template value<E021>() = other.template value<E3>();
    }
    if constexpr (elems::has_scalar(elements))
    {
        out.template value<E0123>() = other.template value<Scalar>();
    }
    return out;
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<first_elements, ScalarType> sandwich_product(const Multivector<first_elements, ScalarType>& object,
                                                         const Multivector<second_elements, ScalarType>& motor)
{
    return Multivector<first_elements>(motor * object * ~motor);
}

template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> conjugate(const Multivector<elements, ScalarType>& object)
{
    using namespace elems;
    Multivector<elements, ScalarType> out;
    if constexpr (has_scalar(elements))
    {
        out.template value<Scalar>() = object.template value<Scalar>();
    }
    if constexpr (has_e0(elements))
    {
        out.template value<E0>() = -object.template value<E0>();
    }
    if constexpr (has_e1(elements))
    {
        out.template value<E1>() = -object.template value<E1>();
    }
    if constexpr (has_e2(elements))
    {
        out.template value<E2>() = -object.template value<E2>();
    }
    if constexpr (has_e3(elements))
    {
        out.template value<E3>() = -object.template value<E3>();
    }
    if constexpr (has_e01(elements))
    {
        out.template value<E01>() = -object.template value<E01>();
    }
    if constexpr (has_e02(elements))
    {
        out.template value<E02>() = -object.template value<E02>();
    }
    if constexpr (has_e03(elements))
    {
        out.template value<E03>() = -object.template value<E03>();
    }
    if constexpr (has_e23(elements))
    {
        out.template value<E23>() = -object.template value<E23>();
    }
    if constexpr (has_e31(elements))
    {
        out.template value<E31>() = -object.template value<E31>();
    }
    if constexpr (has_e12(elements))
    {
        out.template value<E12>() = -object.template value<E12>();
    }
    if constexpr (has_e123(elements))
    {
        out.template value<E123>() = object.template value<E123>();
    }
    if constexpr (has_e032(elements))
    {
        out.template value<E032>() = object.template value<E032>();
    }
    if constexpr (has_e013(elements))
    {
        out.template value<013>() = object.template value<013>();
    }
    if constexpr (has_e021(elements))
    {
        out.template value<E021>() = object.template value<E021>();
    }
    if constexpr (has_e0123(elements))
    {
        out.template value<E0123>() = object.template value<E0123>();
    }
    return out;
};

template <elems::Elems elements, typename ScalarType>
ScalarType norm(const Multivector<elements, ScalarType>& object)
{
    return std::sqrt(std::abs(object * conjugate(object).template values<elems::Scalar>()));
}

/// Pre-defined float-typed primitives
using ScalarF = Multivector<elems::ScalarElems, float>;
using ComplexF = Multivector<elems::ComplexElems, float>;
using PlaneF = Multivector<elems::PlaneElems, float>;
using LineF = Multivector<elems::LineElems, float>;
using PointF = Multivector<elems::PointElems, float>;
using RotorF = Multivector<elems::RotorElems, float>;
using TranslatorF = Multivector<elems::TranslatorElems, float>;
using MotorF = Multivector<elems::MotorElems, float>;

/// Pre-defined double-typed primitives
using ScalarD = Multivector<elems::ScalarElems, double>;
using ComplexD = Multivector<elems::ComplexElems, double>;
using PlaneD = Multivector<elems::PlaneElems, double>;
using LineD = Multivector<elems::LineElems, double>;
using PointD = Multivector<elems::PointElems, double>;
using RotorD = Multivector<elems::RotorElems, double>;
using TranslatorD = Multivector<elems::TranslatorElems, double>;
using MotorD = Multivector<elems::MotorElems, double>;

template <typename ScalarType = float>
static Multivector<elems::RotorElems, ScalarType> Rotor(const ScalarType angle,
                                                        const Multivector<elems::LineElems, ScalarType> line)
{
    return std::cos(angle / ScalarType(2.0)) + std::sin(angle / ScalarType(2.0)) * line.normalized();
}

}  // namespace tiny_pga

#endif  // TINY_PGA_H
