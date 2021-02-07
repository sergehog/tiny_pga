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
#include <gtest/gtest.h>

using namespace tiny_autodf;

TEST(BasicAutodiffTest, OneDependentVariableTest)
{
    AutoDf<> x = 15.F;
    AutoDf<> y = x + 5.F;
    AutoDf<> z = (2.f * x + 2.F) * (y - 3.F);
    AutoDf<> w = z / (x + 1.F);

    EXPECT_EQ(x.value(), 15.f);
    EXPECT_EQ(y.value(), 20.f);
    EXPECT_EQ(z.value(), 544.f);
    EXPECT_EQ(w.value(), 544.f / 16.f);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);
    ASSERT_EQ(z.variables().size(), 1);
    ASSERT_EQ(w.variables().size(), 1);

    y += 5.f;
    EXPECT_EQ(y.value(), 25.f);
    auto xe = x.eval();
    auto ye = y.eval();
    auto ze = z.eval();
    auto we = w.eval();

    EXPECT_EQ(x.value(), xe.value);
    EXPECT_EQ(y.value(), ye.value);
    EXPECT_EQ(z.value(), ze.value);
    EXPECT_EQ(w.value(), we.value);

    ASSERT_EQ(xe.derivatives.size(), 1);
    ASSERT_EQ(ye.derivatives.size(), 1);
    ASSERT_EQ(ze.derivatives.size(), 1);
    ASSERT_EQ(we.derivatives.size(), 1);

    EXPECT_EQ(xe.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ze.derivatives[x.ID()], 66.F);
    EXPECT_EQ(we.derivatives[x.ID()], 2.F);

    // nothing changed, except y
    EXPECT_EQ(x.value(), 15.f);
    EXPECT_EQ(y.value(), 25.f);
    EXPECT_EQ(z.value(), 544.f);
    EXPECT_EQ(w.value(), 544.f / 16.f);

    x.value() += 1.F;
    // re-calculate only one formula, but other are updated too, since they are part of call-graph
    we = w.eval();
    EXPECT_EQ(x.value(), 16.f);
    EXPECT_EQ(y.value(), 25.f);
    EXPECT_EQ(z.value(), 612.f);
    EXPECT_EQ(w.value(), 36.f);
}

TEST(BasicAutodiffTest, SumTest)
{
    AutoDf<float> x = 7.F;
    AutoDf<float> y = (x + 3.F) + 5.F;
    AutoDf<float> z = (5.F + y) + x;

    EXPECT_EQ(x.value(), 7.F);
    EXPECT_EQ(y.value(), 15.F);
    EXPECT_EQ(z.value(), 27.F);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);
    ASSERT_EQ(z.variables().size(), 1);

    x.value() -= 1.F;
    // x value changes, but not others
    EXPECT_EQ(x.value(), 6.F);
    EXPECT_EQ(y.value(), 15.F);
    EXPECT_EQ(z.value(), 27.F);

    // re-evaluate
    auto xe = x.eval();
    auto ye = y.eval();
    auto ze = z.eval();

    // now all values changed
    EXPECT_EQ(x.value(), 6.F);
    EXPECT_EQ(y.value(), 14.F);
    EXPECT_EQ(z.value(), 25.F);

    // evaluation results are the same
    EXPECT_EQ(x.value(), xe.value);
    EXPECT_EQ(y.value(), ye.value);
    EXPECT_EQ(z.value(), ze.value);

    // number of partial derivatives is 1
    ASSERT_EQ(xe.derivatives.size(), 1);
    ASSERT_EQ(ye.derivatives.size(), 1);
    ASSERT_EQ(ze.derivatives.size(), 1);

    // all partial derivatives are 1, since only + was used
    EXPECT_EQ(xe.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ze.derivatives[x.ID()], 2.F);
}

TEST(BasicAutodiffTest, SubtractTest)
{
    AutoDf<> x = 10.F;
    AutoDf<> y = 20.F - x - 5.F;
    AutoDf<> z = 7.F - (10.F - y);

    EXPECT_EQ(x.value(), 10.f);
    EXPECT_EQ(y.value(), 5.F);
    EXPECT_EQ(z.value(), 2.F);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);
    ASSERT_EQ(z.variables().size(), 1);

    x.value() -= 1.F;
    // x value changes, but not others
    EXPECT_EQ(x.value(), 9.f);
    EXPECT_EQ(y.value(), 5.F);
    EXPECT_EQ(z.value(), 2.F);

    // re-evaluate
    auto xe = x.eval();
    auto ye = y.eval();
    auto ze = z.eval();

    // x value changes, but not others
    EXPECT_EQ(x.value(), 9.f);
    EXPECT_EQ(y.value(), 6.F);
    EXPECT_EQ(z.value(), 3.F);

    // evaluation results are the same
    EXPECT_EQ(x.value(), xe.value);
    EXPECT_EQ(y.value(), ye.value);
    EXPECT_EQ(z.value(), ze.value);

    // number of partial derivatives is 1
    ASSERT_EQ(xe.derivatives.size(), 1);
    ASSERT_EQ(ye.derivatives.size(), 1);
    ASSERT_EQ(ze.derivatives.size(), 1);

    // all partial derivatives are 1, since only + was used
    EXPECT_EQ(xe.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID()], -1.F);
    EXPECT_EQ(ze.derivatives[x.ID()], -1.F);
}

TEST(BasicAutodiffTest, MultiplicationTest)
{
    AutoDf<> x = 7.F;
    AutoDf<> y = (x - 1.f) * (x + 1.F) * 2.F;

    EXPECT_EQ(x.value(), 7.f);
    EXPECT_EQ(y.value(), 6.f * 8.F * 2.F);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);

    x.value() -= 1.F;
    // x value changes, but not others
    EXPECT_EQ(x.value(), 6.f);
    EXPECT_EQ(y.value(), 6.f * 8.F * 2.F);

    // re-evaluate
    auto xe = x.eval();
    auto ye = y.eval();

    // x value changes, but not others
    EXPECT_EQ(x.value(), 6.f);
    EXPECT_EQ(y.value(), 5.f * 7.F * 2.F);

    // evaluation results are the same
    EXPECT_EQ(x.value(), xe.value);
    EXPECT_EQ(y.value(), ye.value);

    // number of partial derivatives is 1
    ASSERT_EQ(xe.derivatives.size(), 1);
    ASSERT_EQ(ye.derivatives.size(), 1);

    // all partial derivatives are 1, since only + was used
    EXPECT_EQ(xe.derivatives[x.ID()], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID()], 4 * x.value());
}

TEST(BasicAutodiffTest, DivisionTest)
{
    AutoDf<> x = 7.F;
    AutoDf<> y = (x - 1.f) / (x + 1.F) / 2.F;

    EXPECT_EQ(x.value(), 7.f);
    EXPECT_EQ(y.value(), 6.f / 8.F / 2.F);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);

    x.value() -= 1.F;
    // x value changes, but not others
    EXPECT_EQ(x.value(), 6.f);
    EXPECT_EQ(y.value(), 6.f / 8.F / 2.F);

    // re-evaluate
    auto xe = x.eval();
    auto ye = y.eval();

    // x value changes, but not others
    EXPECT_EQ(x.value(), 6.f);
    EXPECT_EQ(y.value(), 5.f / 7.F / 2.F);

    // evaluation results are the same
    EXPECT_EQ(x.value(), xe.value);
    EXPECT_EQ(y.value(), ye.value);

    // number of partial derivatives is 1
    ASSERT_EQ(xe.derivatives.size(), 1);
    ASSERT_EQ(ye.derivatives.size(), 1);

    // all partial derivatives are 1, since only + was used
    EXPECT_EQ(xe.derivatives[x.ID()], 1.F);
    EXPECT_NEAR(ye.derivatives[x.ID()], 0.0204082F, 0.00001F);
}

TEST(BasicAutodiffTest, AbsMinMaxTest)
{
    AutoDf<>::StartVariables();
    AutoDf<> x = 7.F;
    AutoDf<> y = -5.F;
    AutoDf<>::StartConstants();
    AutoDf<> absx = abs(x);
    AutoDf<> absy = abs(y);

    ASSERT_EQ(absx.value(), 7.F);
    ASSERT_EQ(absy.value(), 5.F);

    ASSERT_EQ(min(x, y).value(), -5.F);
    ASSERT_EQ(max(x, y).value(), 7.F);

    ASSERT_EQ(min(absx, absy).value(), 5.F);
    ASSERT_EQ(max(-absx, -absy).value(), -5.F);

    auto ex = absx.eval();
    auto ey = absy.eval();
    ASSERT_EQ(ex.value, absx.value());
    ASSERT_EQ(ey.value, absy.value());
    ASSERT_EQ(ex.derivatives.size(), 1);
    ASSERT_EQ(ey.derivatives.size(), 1);

    ASSERT_EQ(ex.derivatives.begin()->second, 1.F);
    ASSERT_EQ(ey.derivatives.begin()->second, -1.F);
}

TEST(BasicAutodiffTest, AbsSinCosTest)
{
    AutoDf<>::StartVariables();
    AutoDf<> x = 7.F;
    AutoDf<>::StartConstants();
    AutoDf<> sinx = sin(x);
    AutoDf<> cosx = cos(x);

    ASSERT_EQ(sinx.value(), std::sin(7.F));
    ASSERT_EQ(cosx.value(), std::cos(7.F));

    auto e1 = sinx.eval();
    auto e2 = cosx.eval();

    ASSERT_EQ(e1.value, sinx.value());
    ASSERT_EQ(e2.value, cosx.value());
    ASSERT_EQ(e1.derivatives.size(), 1);
    ASSERT_EQ(e2.derivatives.size(), 1);

    ASSERT_EQ(e1.derivatives.begin()->second, e2.value);
    ASSERT_EQ(e2.derivatives.begin()->second, -e1.value);
}