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

#ifndef TINY_PGA_INNER_PRODUCT_IMPL_H
#define TINY_PGA_INNER_PRODUCT_IMPL_H

/// Inner Product operation
template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::inner_product(first_elements, second_elements), ScalarType> InnerProduct(
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

#endif  // TINY_PGA_INNER_PRODUCT_IMPL_H
