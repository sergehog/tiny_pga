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

#ifndef TINY_PGA_ADDITION_IMPL_H
#define TINY_PGA_ADDITION_IMPL_H

template <elems::Elems first_elements, elems::Elems second_elements, typename ScalarType>
Multivector<elems::addition(first_elements, second_elements), ScalarType> Addition(
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

#endif  // TINY_PGA_ADDITION_IMPL_H
