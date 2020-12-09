//
// Created by Sergey Sergey on 7.12.2020.
//

#include "../tiny_autodf.h"
#include <gtest/gtest.h>

using namespace tiny_autodiff;

TEST(BasicAutodiffTest, OneDependentVariableTest)
{
    AutoDf<> x = 15.F;
    AutoDf<> y = x + 5.F;
    AutoDf<> z = y * 2.F;
    AutoDf<> w = y * z;
    //AutoDf<> h = w / (x+5.F);

    ASSERT_EQ(x.variables().size(), 1);

    EXPECT_EQ(x.value(), 15.f);
    auto xe = x.eval();
    EXPECT_EQ(x.value(), xe.first);
    ASSERT_EQ(xe.second.size(), 1);
    EXPECT_EQ(xe.second[1], 1.F);

    ASSERT_EQ(y.variables().size(), 1);
    EXPECT_EQ(y.value(), 20.f);
    EXPECT_EQ(*(y.variables().begin()->second), 15.f);
    auto ye = y.eval();
    EXPECT_EQ(y.value(), ye.first);
    ASSERT_EQ(ye.second.size(), 1);
    EXPECT_EQ(ye.second[1], 1.F);

    ASSERT_EQ(z.variables().size(), 1);
    EXPECT_EQ(z.value(), 40.f);
    EXPECT_EQ(*(z.variables().begin()->second), 15.f);
    auto ze = z.eval();
    EXPECT_EQ(z.value(), ze.first);
    ASSERT_EQ(ze.second.size(), 1);
    EXPECT_EQ(ze.second[1], 2.F);

    ASSERT_EQ(w.variables().size(), 1);
    EXPECT_EQ(w.value(), 800.f);
    EXPECT_EQ(*(w.variables().begin()->second), 15.f);
    auto we = w.eval();
    EXPECT_EQ(w.value(), we.first);
    ASSERT_EQ(we.second.size(), 1);
    EXPECT_EQ(we.second[1], 80.F);

//    ASSERT_EQ(h.variables().size(), 1);
//    EXPECT_EQ(h.value(), 800.F/15.F);
//    EXPECT_EQ(*(h.variables().begin()->second), x.value());
//    auto he = h.eval();
//    EXPECT_EQ(h.value(), he.first);
//    ASSERT_EQ(he.second.size(), 1);
//    EXPECT_EQ(he.second[1], 400.F);
}
