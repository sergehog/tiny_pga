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
#include "../tiny_pga.h"
#include <gtest/gtest.h>

using namespace tiny_autodf;
using namespace tiny_pga;

using Float = AutoDf<float>;

TEST(AutoDfPGATest, SimpleTest)
{
    Float::StartConstants(false);
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float w = Float(1.F, true);
    Float::StartConstants();

    Multivector<elems::PointElems, Float> X;
    X.e021() = x;
    X.e013() = y;
    X.e032() = z;
    X.e123() = w;

    EXPECT_EQ(X.e021().value(), 2.F);
    EXPECT_EQ(X.e013().value(), 3.F);
    EXPECT_EQ(X.e032().value(), 4.F);
    EXPECT_EQ(X.e123().value(), 1.F);

    Multivector<elems::multiplication(elems::PointElems, elems::PointElems), Float> Xn = X * X;

    EXPECT_EQ(Xn.scalar().value(), -1.F);
    auto variables = Xn.scalar().variables();
    ASSERT_EQ(variables.size(), 0);

    Float formula = (x * y * z + w) / (y - z * x);
    EXPECT_NO_THROW(formula.eval());
    variables = formula.variables();

    ASSERT_EQ(variables.size(), 3);
    EXPECT_EQ(*variables[x.ID()], 2.F);
    EXPECT_EQ(*variables[y.ID()], 3.F);
    EXPECT_EQ(*variables[z.ID()], 4.F);

}