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

#ifndef TINY_PGA_MULTIVECTOR_FUNCTION_DEFINITIONS_H
#define TINY_PGA_MULTIVECTOR_FUNCTION_DEFINITIONS_H

#include "elems.h"

namespace tiny_pga
{

/// Compile-time optimized (using C++ templating) implementation of 3D PGA Multivector
template <elems::Elems elements, typename ScalarType = float>
struct Multivector;

/// Geometric Product (generic implementation)
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::geometric_product(first_elements, second_elements), ScalarType> GeometricProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Inner (dot) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::inner_product(first_elements, second_elements), ScalarType> InnerProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Outer (meet) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::outer_product(first_elements, second_elements), ScalarType> OuterProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Regressive (join) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::regressive_product(first_elements, second_elements), ScalarType> RegressiveProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Commutator (cross) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::commutator_product(first_elements, second_elements), ScalarType> CommutatorProduct(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Sandwich (transform) Product
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<first_elements, ScalarType> SandwichProduct(const Multivector<first_elements, ScalarType>& object,
                                                        const Multivector<second_elements, ScalarType>& transform);

/// Addition
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> Addition(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Subtraction
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> Subtraction(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second);

/// Reverse
template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> Reverse(const Multivector<elements, ScalarType>& object);

/// Dual
template <elems::Elems elements, typename ScalarType>
Multivector<elems::dual(elements), ScalarType> Dual(const Multivector<elements, ScalarType>& object);

/// Clifford Conjugation
template <elems::Elems elements, typename ScalarType>
Multivector<elements, ScalarType> Conjugate(const Multivector<elements, ScalarType>& object);

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::geometric_product(first_elements, second_elements), ScalarType> operator*(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    return tiny_pga::GeometricProduct(first, second);
}

template <elems::Elems first_elements, typename ScalarType>
inline Multivector<elems::geometric_product(first_elements, elems::Elems(elems::Values::kScalar)), ScalarType>
operator*(const Multivector<first_elements, ScalarType>& first, const ScalarType& scalar)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> second{scalar};
    return tiny_pga::GeometricProduct(first, second);
}

template <elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::geometric_product(elems::Elems(elems::Values::kScalar), second_elements), ScalarType>
operator*(const ScalarType& scalar, const Multivector<second_elements, ScalarType>& second)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> first{scalar};
    return tiny_pga::GeometricProduct(first, second);
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, second_elements), ScalarType> operator+(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    return tiny_pga::Addition(first, second);
}

template <elems::Elems first_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, elems::Elems(elems::Values::kScalar)), ScalarType> operator+(
    const Multivector<first_elements, ScalarType>& first,
    const ScalarType& scalar)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> second{scalar};
    return tiny_pga::Addition(first, second);
}

template <elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(elems::Elems(elems::Values::kScalar), second_elements), ScalarType> operator+(
    const ScalarType& scalar,
    const Multivector<second_elements, ScalarType>& second)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> first{scalar};
    return tiny_pga::Addition(first, second);
}

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, second_elements), ScalarType> operator-(
    const Multivector<first_elements, ScalarType>& first,
    const Multivector<second_elements, ScalarType>& second)
{
    return tiny_pga::Subtraction(first, second);
}

template <elems::Elems second_elements, typename ScalarType>
inline Multivector<elems::addition(elems::Elems(elems::Values::kScalar), second_elements), ScalarType> operator-(
    const ScalarType& scalar,
    const Multivector<second_elements, ScalarType>& second)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> first{scalar};
    return tiny_pga::Subtraction(first, second);
}

template <elems::Elems first_elements, typename ScalarType>
inline Multivector<elems::addition(first_elements, elems::Elems(elems::Values::kScalar)), ScalarType> operator-(
    const Multivector<first_elements, ScalarType>& first,
    const ScalarType& scalar)
{
    Multivector<elems::Elems(elems::Values::kScalar), ScalarType> second{scalar};
    return tiny_pga::Subtraction(first, second);
}
} // namespace tiny_pga

#endif  // TINY_PGA_MULTIVECTOR_FUNCTION_DEFINITIONS_H
