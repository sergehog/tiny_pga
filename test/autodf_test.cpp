//
// Created by Sergey Sergey on 7.12.2020.
//

#include "../tiny_autodf.h"
#include <gtest/gtest.h>

using namespace tiny_autodf;

TEST(BasicAutodiffTest, OneDependentVariableTest)
{
    AutoDf<> x = 15.F;
    AutoDf<> y = x + 5.F;
    AutoDf<> z = (2.f * x + 2.F) * (y - 3.F);
    AutoDf<> w = 0.5 * z / (x + 1.F);

    EXPECT_EQ(x.value(), 15.f);
    EXPECT_EQ(y.value(), 20.f);
    EXPECT_EQ(z.value(), 544.f);
    EXPECT_EQ(w.value(), 0.5f * 544.f / 16.f);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);
    ASSERT_EQ(z.variables().size(), 1);
    ASSERT_EQ(w.variables().size(), 1);

    auto xe = x.eval();
    auto ye = y.eval();
    auto ze = z.eval();
    auto we = w.eval();

    EXPECT_EQ(x.value(), xe.value);
    ASSERT_EQ(xe.derivatives.size(), 1);
    EXPECT_EQ(xe.derivatives[x.ID], 1.F);

    EXPECT_EQ(y.value(), ye.value);
    ASSERT_EQ(ye.derivatives.size(), 1);
    EXPECT_EQ(ye.derivatives[x.ID], 1.F);

    EXPECT_EQ(z.value(), ze.value);
    ASSERT_EQ(ze.derivatives.size(), 1);
    EXPECT_EQ(ze.derivatives[x.ID], 66.F);

    EXPECT_EQ(w.value(), we.value);
    ASSERT_EQ(we.derivatives.size(), 1);
    EXPECT_EQ(we.derivatives[x.ID], 256.F);
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
    EXPECT_EQ(xe.derivatives[x.ID], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID], 1.F);
    EXPECT_EQ(ze.derivatives[x.ID], 2.F);
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
    EXPECT_EQ(xe.derivatives[x.ID], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID], -1.F);
    EXPECT_EQ(ze.derivatives[x.ID], -1.F);
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
    EXPECT_EQ(xe.derivatives[x.ID], 1.F);
    EXPECT_EQ(ye.derivatives[x.ID], 4 * x.value());
}