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


#include "../tiny_autodf.h"
#include "../pga3d.h"
#include <gtest/gtest.h>

using namespace tiny_autodf;

using Float = AutoDf<float>;
using APGA = PGA3D<Float>;
using PGA = PGA3D<float>;

TEST(AutoDfPGA3DTest, SimpleTest)
{
    Float::SetType(Float::AutoType::kVariableType);
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float::SetType(Float::AutoType::kConstType);

    const APGA e0(kE0);
    APGA aa = e0 * x;

    APGA p = point(x, y, z);
//
    auto p_sq = p*p;
    auto vars1 = p_sq[0].variables();
    EXPECT_EQ(vars1.size(), 3);
    EXPECT_EQ(p_sq[0].value(), -1);
//    x = 3.F;
//    auto vars2 = p_sq[0].eval();
//    EXPECT_EQ(vars2.derivatives.size(), 3);


//    Multivector<elems::multiplication(elems::PointElems, elems::PointElems), Float> Xn = X * X;
//    EXPECT_EQ(Xn.scalar(), -1.F);
//    EXPECT_EQ(Xn.scalar().value(), -1.F);
//    auto variables = Xn.scalar().variables();
//    EXPECT_EQ(variables.size(), 3);
//    auto xv = variables[x.ID];
//    std::cout << *xv << std::endl;

//    EXPECT_EQ(variables[x.ID], x.value());
//    EXPECT_EQ(variables[y.ID], y.value());
//    EXPECT_EQ(variables[z.ID], z.value());


}