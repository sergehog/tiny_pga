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
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float w = Float(1.F, true);

    Multivector<elems::PointElems, Float> X {{},{},{}, {x, y, z, w}};
    EXPECT_EQ(X.e123().value(), 1.F);
    EXPECT_EQ(X.e021().value(), 2.F);
    EXPECT_EQ(X.e013().value(), 3.F);
    EXPECT_EQ(X.e032().value(), 4.F);

    Multivector<elems::multiplication(elems::PointElems, elems::PointElems), Float> Xn = X * X;
    EXPECT_EQ(Xn.scalar(), -1.F);
    EXPECT_EQ(Xn.scalar().value(), -1.F);
    auto variables = Xn.scalar().variables();
    EXPECT_EQ(variables.size(), 3);
    auto xv = variables[x.ID];
    std::cout << *xv << std::endl;

//    EXPECT_EQ(variables[x.ID], x.value());
//    EXPECT_EQ(variables[y.ID], y.value());
//    EXPECT_EQ(variables[z.ID], z.value());


}