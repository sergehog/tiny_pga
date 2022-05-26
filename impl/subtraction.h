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
#ifndef TINY_PGA_SUBTRACTION_IMPL_H
#define TINY_PGA_SUBTRACTION_IMPL_H

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> Subtraction(
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

#endif  // TINY_PGA_SUBTRACTION_IMPL_H
