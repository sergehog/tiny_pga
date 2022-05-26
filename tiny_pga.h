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

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <type_traits>

#include "include/multivector_function_definitions.h"

namespace tiny_pga
{

template <elems::Elems elements, typename ScalarType>
struct Multivector
{
    /// Exposing elements outside, as static const value
    static constexpr elems::Elems Elements = elements;

    /// values of Multivector elements
    std::array<ScalarType, elems::count(elements)> values{};

    /// "Shortified" Ctor - USE WITH CAUTION! Length and order of elements directly depend of type of Multivector
    Multivector(std::initializer_list<ScalarType> list) { std::copy(list.begin(), list.end(), values.begin()); }

    /// explicit cast from other Multivector type
    template <elems::Elems other_elements, typename OtherScalarType>
    explicit Multivector(const Multivector<other_elements, OtherScalarType>& other)
    {
#define ELEM_CAST(elem)                                                                     \
    if constexpr (elems::has_elem<elem>(elements) && elems::has_elem<elem>(other_elements)) \
    {                                                                                       \
        value<elem>() = static_cast<ScalarType>(other.template value<elem>());              \
    }
        ELEM_CAST(elems::Scalar);
        ELEM_CAST(elems::E0);
        ELEM_CAST(elems::E1);
        ELEM_CAST(elems::E2);
        ELEM_CAST(elems::E3);
        ELEM_CAST(elems::E01);
        ELEM_CAST(elems::E02);
        ELEM_CAST(elems::E03);
        ELEM_CAST(elems::E12);
        ELEM_CAST(elems::E31);
        ELEM_CAST(elems::E23);
        ELEM_CAST(elems::E021);
        ELEM_CAST(elems::E013);
        ELEM_CAST(elems::E032);
        ELEM_CAST(elems::E123);
        ELEM_CAST(elems::E0123);
#undef ELEM_CAST
    }

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

    /// Geometric Product operator
    template <elems::Elems other_elements>
    Multivector<elems::geometric_product(elements, other_elements), ScalarType> operator*(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::GeometricProduct(*this, other);
    }

    /// Inner product operator (also known as Dot-product)
    template <elems::Elems other_elements>
    Multivector<elems::inner_product(elements, other_elements), ScalarType> operator|(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::InnerProduct(*this, other);
    }

    /// Alias for Inner product (Dot-product)
    template <elems::Elems other_elements>
    Multivector<elems::inner_product(elements, other_elements), ScalarType> Dot(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::InnerProduct(*this, other);
    }

    /// Outer product operator (also known as Wedge- or Meet- product)
    template <elems::Elems other_elements>
    Multivector<elems::outer_product(elements, other_elements), ScalarType> operator^(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::OuterProduct(*this, other);
    }

    /// Alias for Outer product operator
    template <elems::Elems other_elements>
    Multivector<elems::outer_product(elements, other_elements), ScalarType> Meet(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::OuterProduct(*this, other);
    }

    /// Regressive product operator (also known as Vee- or Join-product)
    template <elems::Elems other_elements>
    Multivector<elems::regressive_product(elements, other_elements), ScalarType> operator&(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::RegressiveProduct(*this, other);
    }

    /// Alias for Regressive product operator
    template <elems::Elems other_elements>
    Multivector<elems::regressive_product(elements, other_elements), ScalarType> Join(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::RegressiveProduct(*this, other);
    }

    /// Multivector addition
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> operator+(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::Addition(*this, other);
    }

    /// Alias for addition
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> Add(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::Addition(*this, other);
    }

    /// Multivector subtraction
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> operator-(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::Subtraction(*this, other);
    }

    /// Unary minus operator
    Multivector<elements, ScalarType> operator-() const
    {
        Multivector<0, ScalarType> nothing{};
        return tiny_pga::Subtraction(nothing, *this);
    }

    /// Alias for subtraction
    template <elems::Elems other_elements>
    Multivector<elems::addition(elements, other_elements), ScalarType> Sub(
        const Multivector<other_elements, ScalarType>& other) const
    {
        return tiny_pga::Subtraction(*this, other);
    }

    /// Reverse operator
    Multivector<elements, ScalarType> operator~() const
    {
        return tiny_pga::Reverse(*this);
    }

    /// Reverse function
    Multivector<elements, ScalarType> Reverse() const
    {
        return tiny_pga::Reverse(*this);
    }

    /// Duality operator
    Multivector<elems::dual(elements), ScalarType> operator!() const
    {
        return tiny_pga::Dual(*this);
    }

    /// Duality function
    Multivector<elems::dual(elements), ScalarType> Dual() const
    {
        return tiny_pga::Dual(*this);
    }

    /// sandwich product with Rotor
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::RotorElems, ScalarType>& other) const
    {
        return tiny_pga::SandwichProduct(*this, other);
    }

    /// sandwich product with Translator
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::TranslatorElems, ScalarType>& other) const
    {
        return tiny_pga::SandwichProduct(*this, other);
    }

    /// sandwich product with Motor
    Multivector<elements, ScalarType> operator<<(const Multivector<elems::MotorElems, ScalarType>& other) const
    {
        return tiny_pga::SandwichProduct(*this, other);
    }

    Multivector<elements, ScalarType> Normalized() const
    {
        return *this * ScalarType(1.F / Norm(*this));
    }
};

#include "impl/geometric_propduct.hpp"
#include "impl/inner_product.hpp"
#include "impl/outer_product.hpp"
#include "impl/addition.hpp"
#include "impl/subtraction.h"


template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> Reverse(const Multivector<elements, ScalarType>& object)
{
    using namespace elems;
    Multivector<elements, ScalarType> out{};
    if constexpr (has_scalar(elements))
    {
        out.template value<Scalar>() = object.template value<Scalar>();
    }
    if constexpr (has_e0(elements))
    {
        out.template value<E0>() = object.template value<E0>();
    }
    if constexpr (has_e1(elements))
    {
        out.template value<E1>() = object.template value<E1>();
    }
    if constexpr (has_e2(elements))
    {
        out.template value<E2>() = object.template value<E2>();
    }
    if constexpr (has_e3(elements))
    {
        out.template value<E3>() = object.template value<E3>();
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
    if constexpr (has_e12(elements))
    {
        out.template value<E12>() = -object.template value<E12>();
    }
    if constexpr (has_e31(elements))
    {
        out.template value<E31>() = -object.template value<E31>();
    }
    if constexpr (has_e23(elements))
    {
        out.template value<E23>() = -object.template value<E23>();
    }
    if constexpr (has_e021(elements))
    {
        out.template value<E021>() = -object.template value<E021>();
    }
    if constexpr (has_e013(elements))
    {
        out.template value<E013>() = -object.template value<E013>();
    }
    if constexpr (has_e032(elements))
    {
        out.template value<E032>() = -object.template value<E032>();
    }
    if constexpr (has_e123(elements))
    {
        out.template value<E123>() = -object.template value<E123>();
    }
    if constexpr (has_e0123(elements))
    {
        out.template value<E0123>() = object.template value<E0123>();
    }
    return out;
}

template <elems::Elems elements, typename ScalarType>
Multivector<elems::dual(elements), ScalarType> Dual(const Multivector<elements, ScalarType>& other)
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
Multivector<first_elements, ScalarType> SandwichProduct(const Multivector<first_elements, ScalarType>& object,
                                                         const Multivector<second_elements, ScalarType>& transform)
{
    return Multivector<first_elements>(transform * object * ~transform);
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::regressive_product(first_elements, second_elements), ScalarType> RegressiveProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    constexpr auto type = elems::regressive_product(first_elements, second_elements);
    return Multivector<type, ScalarType>(Dual(Dual(first) ^ Dual(second)));
}

template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> Conjugate(const Multivector<elements, ScalarType>& object)
{
    using namespace elems;
    Multivector<elements, ScalarType> out{};
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
ScalarType Norm(const Multivector<elements, ScalarType>& object)
{
    return std::sqrt(std::abs((object * Conjugate(object)).template value<elems::Scalar>()));
}

/// Rotor around Euclidean line
template <typename ScalarType = float>
static Multivector<elems::RotorElems, ScalarType> Rotor(const ScalarType angle,
                                                        const Multivector<elems::LineElems, ScalarType>& line)
{
    return Multivector<elems::RotorElems, ScalarType>(std::cos(angle / ScalarType(2.0)) +
                                                      std::sin(angle / ScalarType(2.0)) * line.Normalized());
}

/// Translator over directed line
template <typename ScalarType = float>
static Multivector<elems::TranslatorElems, ScalarType> Translator(const ScalarType dist,
                                                                  const Multivector<elems::LineElems, ScalarType>& line)
{
    return ScalarType(1) + (dist / ScalarType(2.0)) * line;
}

/// Translator
template <typename ScalarType = float>
static Multivector<elems::TranslatorElems, ScalarType> Translator(const ScalarType dx,
                                                                  const ScalarType dy,
                                                                  const ScalarType dz)
{
    return Multivector<elems::PlaneElems, ScalarType>{1, -dx / 2., -dy / 2., -dz / 2.};
}

/// A plane is defined using its homogenous equation ax + by + cz + d = 0
template <typename ScalarType = float>
static Multivector<elems::PlaneElems, ScalarType> Plane(ScalarType a, ScalarType b, ScalarType c, ScalarType d)
{
    return Multivector<elems::PlaneElems, ScalarType>{d, a, b, c};
}

/// A point is just a homogeneous point, euclidean coordinates plus the origin
template <typename ScalarType = float>
static Multivector<elems::PointElems, ScalarType> Point(const ScalarType& x, const ScalarType& y, const ScalarType& z)
{
    return Multivector<elems::PointElems, ScalarType>{x, y, z, 1.};
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

namespace float_basis
{
const static Multivector<elems::Elems(elems::Values::kE0), float> e0{1};
const static Multivector<elems::Elems(elems::Values::kE1), float> e1{1};
const static Multivector<elems::Elems(elems::Values::kE2), float> e2{1};
const static Multivector<elems::Elems(elems::Values::kE3), float> e3{1};
const static Multivector<elems::Elems(elems::Values::kE01), float> e01{1};
const static Multivector<elems::Elems(elems::Values::kE02), float> e02{1};
const static Multivector<elems::Elems(elems::Values::kE03), float> e03{1};
const static Multivector<elems::Elems(elems::Values::kE12), float> e12{1};
const static Multivector<elems::Elems(elems::Values::kE23), float> e23{1};
const static Multivector<elems::Elems(elems::Values::kE31), float> e31{1};
const static Multivector<elems::Elems(elems::Values::kE021), float> e021{1};
const static Multivector<elems::Elems(elems::Values::kE032), float> e032{1};
const static Multivector<elems::Elems(elems::Values::kE013), float> e013{1};
const static Multivector<elems::Elems(elems::Values::kE123), float> e123{1};
const static Multivector<elems::Elems(elems::Values::kE0123), float> e0123{1};
}  // namespace float_basis
}  // namespace tiny_pga

#endif  // TINY_PGA_H
