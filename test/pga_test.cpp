/*
 * This file is part of the Tiny-PGA distribution (https://github.com/sergehog/tiny_pga)
 * Copyright (c) 2020-2021 Sergey Smirnov / Seregium Oy.
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

using namespace tiny_pga;
using namespace tiny_pga::elems;

class PgaTest : public ::testing::Test
{
  public:
    PgaTest() = default;
};

/// Multivector works as a scalar
TEST_F(PgaTest, ScalarTest)
{
    ScalarF mv{11.F};
    EXPECT_TRUE(elems::has_scalar(mv.Elements));
    EXPECT_EQ(mv[Scalar], 11.F);
    EXPECT_EQ(mv.value<Scalar>(), 11.F);

    mv[Scalar] = 15.f;
    EXPECT_EQ(mv[Scalar], 15.F);
    EXPECT_EQ(mv.value<Scalar>(), 15.F);

    auto mv2 = mv * mv;
    EXPECT_EQ(mv2.Elements, mv.Elements);
    EXPECT_TRUE(elems::has_scalar(mv2.Elements));
    EXPECT_EQ(mv2[Scalar], 15.F * 15.F);
    EXPECT_EQ(mv2.value<Scalar>(), 15.F * 15.F);
}

/// Multivector works as a complex number
TEST_F(PgaTest, ComplexTest)
{
    ComplexF mv{0.F, 1.F};
    EXPECT_TRUE(elems::has_scalar(mv.Elements));
    EXPECT_TRUE(elems::has_e12(mv.Elements));

    EXPECT_EQ(mv[Scalar], 0.F);
    EXPECT_EQ(mv.value<Scalar>(), 0.F);
    EXPECT_EQ(mv[E12], 1.F);
    EXPECT_EQ(mv.value<E12>(), 1.F);

    auto mv2 = mv * mv;
    EXPECT_EQ(mv2.Elements, mv.Elements);
    EXPECT_TRUE(elems::has_scalar(mv2.Elements));
    EXPECT_EQ(mv2[Scalar], -1.F);
    EXPECT_EQ(mv2.value<Scalar>(), -1.F);
    EXPECT_EQ(mv2[E12], 0.F);
    EXPECT_EQ(mv2.value<E12>(), 0.F);
}

TEST_F(PgaTest, RotorTest)
{
    RotorF rotor{1.F, 0.F, 0.F, 0.F};
    PointF point{1.f, 2.f, 3.f, 1.F};

    // Rotor is identity, so point2 values must be the same as in point
    const auto point2 = rotor * point * ~rotor;
    const auto point3 = point << rotor;

    EXPECT_NEAR(point2[E021], point[E021], 1e-5);
    EXPECT_NEAR(point2[E013], point[E013], 1e-5);
    EXPECT_NEAR(point2[E032], point[E032], 1e-5);
    EXPECT_NEAR(point2[E123], point[E123], 1e-5);

    EXPECT_EQ(point3.Elements, point.Elements);
}

TEST_F(PgaTest, PlaneTest)
{
    // Plane is pure vector in PGA
    PlaneF p{1, 2, 3, 4};
    // squaring with Dot-product
    auto p_dot = p | p;  // must be a scalar
    EXPECT_EQ(p_dot.Elements, ScalarElems);

    // squaring with geometric product
    auto p_mul = p * p;
    EXPECT_GT(p_mul.Elements, ScalarElems);  // scalar + all vectors + all bi-vectors

    // squaring both ways must result in same scalar
    EXPECT_EQ(p_mul[Scalar], p_dot[Scalar]);

    auto p_0 = p_mul - p_dot;
    EXPECT_EQ(p_mul.Elements, p_0.Elements);
    EXPECT_EQ(p_0[Scalar], 0.);

    auto p_dual = !p_0;
    EXPECT_NE(p_dual.Elements, p_0.Elements);
    EXPECT_TRUE(elems::has_e0123(p_dual.Elements));
    EXPECT_EQ(p_dual[E0123], 0.);

    auto p_all = p_0 + p_dual;
}

TEST(BasicTest, RotorTest)
{
    RotorF R{0.8F, .2F, .3F, .4F};
    PointF p{1.f, 2.f, 3.f, 4.F};

    // rotated point
    auto p2 = R * p * ~R;

    // rotated back
    auto p3 = ~R * p2 * R;

    //    EXPECT_NEAR(p[E021], p3[E021], 1e-5);
    //    EXPECT_NEAR(p[E013], p3[E013], 1e-5);
    //    EXPECT_NEAR(p[E032], p3[E032], 1e-5);
    //    EXPECT_NEAR(p[E123], p3[E123], 1e-5);
}

// Check if geometric product rule (A*B = A|B + A^B) holds for vectors
TEST_F(PgaTest, GeometricProductOfVectorIsSumOfInnerAndOuterProductsTest)
{
    PlaneF A{1, 2, 3, 4};
    PlaneF B{4, 3, 2, 1};

    auto GP = A * B;
    auto IP = A | B;
    auto OP = A ^ B;
    auto GP2 = IP + OP;
    auto Diff = GP - GP2;

    EXPECT_EQ(GP.Elements, GP2.Elements);
    for (auto v : Diff.values)
    {
        EXPECT_NEAR(v, 0.F, 1e-8);
    }
}

/// Test class for checking all possible Multivector variants
template <typename T>
class AllMultivectorsTest : public ::testing::Test
{
  public:
    AllMultivectorsTest() = default;
};

TYPED_TEST_SUITE_P(AllMultivectorsTest);

TYPED_TEST_P(AllMultivectorsTest, BasicTest)
{
    // Inside a test, refer to TypeParam to get the type parameter.
    TypeParam mv{};
    for (std::size_t i = 0UL; i < mv.values.size(); i++)
    {
        mv.values[i] = i;
    }

    auto gp = mv * mv;  // geometric product
    auto ip = mv | mv;  // inner product
    auto op = mv ^ mv;  // wedge product

    EXPECT_GE(gp.Elements, ip.Elements);
    EXPECT_GE(gp.Elements, op.Elements);
    EXPECT_EQ(elems::has_scalar(gp.Elements), elems::has_scalar(ip.Elements));

    if (elems::has_scalar(gp.Elements))
    {
        EXPECT_EQ(gp[Scalar], ip[Scalar]);
    }
}

REGISTER_TYPED_TEST_SUITE_P(AllMultivectorsTest, BasicTest);

using MultivectorTypes =
    ::testing::Types<Multivector<Scalar>, PlaneF, ComplexF, LineF, PointF, RotorF, TranslatorF, MotorF>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, AllMultivectorsTest, MultivectorTypes);
