//
// Copyright Sergey Smirnov  / Seregium Oy 2020
// Email: sergei.smirnov@gmail.com
//

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