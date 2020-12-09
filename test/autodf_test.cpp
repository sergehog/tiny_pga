//
// Created by Sergey Sergey on 7.12.2020.
//

#include "../tiny_autodf.h"
#include <gtest/gtest.h>

using namespace tiny_autodiff;
//
//TEST(BasicAutodiffTest, OneDependentVariableTest)
//{
//    AutoDf<> x = 15.F;
//    AutoDf<> y = x + 5.F;
//    AutoDf<> z = y * 2.F;
//    AutoDf<> w = y * z;
//    //AutoDf<> h = w / (x+5.F);
//
//    ASSERT_EQ(x.variables().size(), 1);
//
//    EXPECT_EQ(x.value(), 15.f);
//    auto xe = x.eval();
//    EXPECT_EQ(x.value(), xe.value);
//    ASSERT_EQ(xe.derivatives.size(), 1);
//    EXPECT_EQ(xe.derivatives[1U], 1.F);
//
//    ASSERT_EQ(y.variables().size(), 1);
//    EXPECT_EQ(y.value(), 20.f);
//    EXPECT_EQ(*(y.variables().begin()->second), 15.f);
//    auto ye = y.eval();
//    EXPECT_EQ(y.value(), ye.value);
//    ASSERT_EQ(ye.derivatives.size(), 1);
//    EXPECT_EQ(ye.derivatives[1U], 1.F);
//
//    ASSERT_EQ(z.variables().size(), 1);
//    EXPECT_EQ(z.value(), 40.f);
//    EXPECT_EQ(*(z.variables().begin()->second), 15.f);
//    auto ze = z.eval();
//    EXPECT_EQ(z.value(), ze.value);
//    ASSERT_EQ(ze.derivatives.size(), 1);
//    EXPECT_EQ(ze.derivatives[1U], 2.F);
//
//    ASSERT_EQ(w.variables().size(), 1);
//    EXPECT_EQ(w.value(), 800.f);
//    EXPECT_EQ(*(w.variables().begin()->second), 15.f);
//    auto we = w.eval();
//    EXPECT_EQ(w.value(), we.value);
//    ASSERT_EQ(we.derivatives.size(), 1);
//    EXPECT_EQ(we.derivatives[1], 80.F);
//
////    ASSERT_EQ(h.variables().size(), 1);
////    EXPECT_EQ(h.value(), 800.F/15.F);
////    EXPECT_EQ(*(h.variables().begin()->second), x.value());
////    auto he = h.eval();
////    EXPECT_EQ(h.value(), he.first);
////    ASSERT_EQ(he.second.size(), 1);
////    EXPECT_EQ(he.second[1], 400.F);
//}


TEST(BasicAutodiffTest, SumTest)
{
    constexpr float kX = 7.F;
    constexpr float kY = 11.F;

    AutoDf<float> x = kX;
    AutoDf<float> y = x + kY;
    AutoDf<float> z = (5.F + y) + x;

    EXPECT_EQ(x.value(), kX);
    EXPECT_EQ(y.value(), kX+kY);
    EXPECT_EQ(z.value(), 5.F+2*kX+kY);

    ASSERT_EQ(x.variables().size(), 1);
    ASSERT_EQ(y.variables().size(), 1);
    ASSERT_EQ(z.variables().size(), 1);

    x.value() -= 1.F;
    // x value changes, but not others
    EXPECT_EQ(x.value(), kX-1);
    EXPECT_EQ(y.value(), kX+kY);
    EXPECT_EQ(z.value(), 5.F+2*kX+kY);

    // re-evaluate
    auto xe = x.eval();
    auto ye = y.eval();
    auto ze = z.eval();

    // now all values changed
    EXPECT_EQ(x.value(), kX-1);
    EXPECT_EQ(y.value(), kX-1+kY);
    EXPECT_EQ(z.value(), 5.F + 2*(kX-1.F) + kY);

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
//
//TEST(BasicAutodiffTest, SubtractTest)
//{
//    AutoDf<> x = 7.F;
//    AutoDf<> y = x - 11.F;
//    AutoDf<> z = 25.F - y;
//
//    EXPECT_EQ(x.value(), 7.f);
//    EXPECT_EQ(y.value(), 7.f - 11.F);
//    EXPECT_EQ(z.value(), 25.f - 7.f + 11.F);
//
//    ASSERT_EQ(x.variables().size(), 1);
//    ASSERT_EQ(y.variables().size(), 1);
//    ASSERT_EQ(z.variables().size(), 1);
//
//    x.value() -= 1.F;
//    // x value changes, but not others
//    EXPECT_EQ(x.value(), 6.f);
//    EXPECT_EQ(y.value(), 7.f - 11.F);
//    EXPECT_EQ(z.value(), 25.f - 7.f + 11.F);
//
//    // re-evaluate
//    auto xe = x.eval();
//    auto ye = y.eval();
//    auto ze = z.eval();
//
//    // x value changes, but not others
//    EXPECT_EQ(x.value(), 6.f);
//    EXPECT_EQ(y.value(), 6.f - 11.F);
//    EXPECT_EQ(z.value(), 25.f - 6.f + 11.F);
//
//    // evaluation results are the same
//    EXPECT_EQ(x.value(), xe.value);
//    EXPECT_EQ(y.value(), ye.value);
//    EXPECT_EQ(z.value(), ze.value);
//
//    // number of partial derivatives is 1
//    ASSERT_EQ(xe.derivatives.size(), 1);
//    ASSERT_EQ(ye.derivatives.size(), 1);
//    ASSERT_EQ(ze.derivatives.size(), 1);
//
//    // all partial derivatives are 1, since only + was used
//    EXPECT_EQ(xe.derivatives[x.ID], 1.F);
//    EXPECT_EQ(ye.derivatives[x.ID], 1.F);
//    EXPECT_EQ(ze.derivatives[x.ID], -1.F);
//}
//
//TEST(BasicAutodiffTest, MultiplicationTest)
//{
//    AutoDf<> x = 7.F;
//    AutoDf<> y = (x - 1.f) * (x + 1.F) * 2.F;
//
//    EXPECT_EQ(x.value(), 7.f);
//    EXPECT_EQ(y.value(), 6.f*8.F*2.F);
//
//    ASSERT_EQ(x.variables().size(), 1);
//    ASSERT_EQ(y.variables().size(), 1);
//
//    x.value() -= 1.F;
//    // x value changes, but not others
//    EXPECT_EQ(x.value(), 6.f);
//    EXPECT_EQ(y.value(), 6.f*8.F*2.F);
//
//    // re-evaluate
//    auto xe = x.eval();
//    auto ye = y.eval();
//
//    // x value changes, but not others
//    EXPECT_EQ(x.value(), 6.f);
//    EXPECT_EQ(y.value(), 5.f*7.F*2.F);
//
//    // evaluation results are the same
//    EXPECT_EQ(x.value(), xe.value);
//    EXPECT_EQ(y.value(), ye.value);
//
//    // number of partial derivatives is 1
//    ASSERT_EQ(xe.derivatives.size(), 1);
//    ASSERT_EQ(ye.derivatives.size(), 1);
//
//
//    // all partial derivatives are 1, since only + was used
//    EXPECT_EQ(xe.derivatives[x.ID], 1.F);
//    EXPECT_EQ(ye.derivatives[x.ID], 4*x.value());
//}