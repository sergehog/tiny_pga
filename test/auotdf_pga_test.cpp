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

#include "../tiny_pga.h"
#include <gtest/gtest.h>
#include <tiny_autodf/tiny_autodf.h>

using namespace tiny_autodf;
using namespace tiny_pga;
using namespace tiny_pga::elems;

using Float = AutoDf<float>;

TEST(AutoDfPGATest, SimpleTest)
{
    Float::VariablesByDefault();
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float w = Float(1.F, true);
    Float::ConstantsByDefault();

    Multivector<elems::PointElems, Float> X{};
    X[E021] = x;
    X[E013] = y;
    X[E032] = z;
    X[E123] = w;

    EXPECT_EQ(X[E021].value(), 2.F);
    EXPECT_EQ(X[E013].value(), 3.F);
    EXPECT_EQ(X[E032].value(), 4.F);
    EXPECT_EQ(X[E123].value(), 1.F);

    auto Xn = X * X;

    EXPECT_EQ(Xn[Scalar].value(), -1.F);
    auto variables = Xn[Scalar].variables();
    ASSERT_EQ(variables.size(), 0);

    Float formula = (x * y * z + w) / (y - z * x);
    EXPECT_NO_THROW(formula.eval());
    variables = formula.variables();

    ASSERT_EQ(variables.size(), 3);
    EXPECT_EQ(*variables[x.ID()], 2.F);
    EXPECT_EQ(*variables[y.ID()], 3.F);
    EXPECT_EQ(*variables[z.ID()], 4.F);
}