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

TEST_F(PgaTest, BasicElemTest)
{
    constexpr Elems elem_1 = static_cast<Elems>(Values::kScalar);
    EXPECT_TRUE(elems::has_scalar(elem_1));

    Multivector<elem_1> mv{11.F};
    EXPECT_TRUE(elems::has_scalar(mv.Elements));
    EXPECT_EQ(mv[Scalar], 11.F);
    EXPECT_EQ(mv.value<Scalar>(), 11.F);

    mv[Scalar] = 15.f;
    EXPECT_EQ(mv[Scalar], 15.F);
    EXPECT_EQ(mv.value<Scalar>(), 15.F);

    auto mv2 = mv * mv;
    EXPECT_EQ(mv2.Elements, mv.Elements);
    EXPECT_TRUE(elems::has_scalar(mv2.Elements));
    EXPECT_EQ(mv2[Scalar], 15.F*15.F);
    EXPECT_EQ(mv2.value<Scalar>(), 15.F*15.F);

    // Plane is pure vector in PGA
    PlaneF p{1, 2, 3, 4};
    // squaring with dot-product
    auto p_dot = p & p; // must be a scalar
    EXPECT_EQ(p_dot.Elements, elem_1);

    // squaring with geometric product
    auto p_mul = p * p;
    EXPECT_GT(p_mul.Elements, elem_1); // scalar + all vectors + all bi-vectors

    // squaring both ways must result in same scalar
    EXPECT_EQ(p_mul[Scalar], p_dot[Scalar]);

    auto p_0 = p_mul - p_dot;
    EXPECT_EQ(p_mul.Elements, p_0.Elements);
    EXPECT_EQ(p_0[Scalar],0.);

    auto p_dual = !p_0;
    EXPECT_NE(p_dual.Elements, p_0.Elements);
    EXPECT_TRUE(elems::has_e0123(p_dual.Elements));
    EXPECT_EQ(p_dual[E0123],0.);

    auto p_all = p_0 + p_dual;
}

/// Test class for checking all possible Multivector variants
template<typename T>
class AllMultivectorsTest
    : public ::testing::Test
{
  public:
    AllMultivectorsTest() = default;
};

TYPED_TEST_SUITE_P(AllMultivectorsTest);

TYPED_TEST_P(AllMultivectorsTest, BasicTest)
{
    // Inside a test, refer to TypeParam to get the type parameter.
    TypeParam mv{};
    for(std::size_t i=0UL; i<mv.values.size(); i++)
    {
        mv.values[i] = i;
    }

    auto mv2 = mv * mv;
    auto mv3 = mv & mv;


    EXPECT_GE(mv2.Elements, mv3.Elements);
    EXPECT_EQ(elems::has_scalar(mv2.Elements), elems::has_scalar(mv3.Elements));

    if(elems::has_scalar(mv2.Elements))
    {
        EXPECT_EQ(mv2[Scalar], mv3[Scalar]);
    }
}

REGISTER_TYPED_TEST_SUITE_P(AllMultivectorsTest,
                            BasicTest);

using MultivectorTypes = ::testing::Types<Multivector<Scalar>, PlaneF, ComplexF, LineF, PointF, RotorF, TranslatorF, MotorF>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, AllMultivectorsTest, MultivectorTypes);



//
//TEST(BasicTest, SquaringTest)
//{
//    tiny_pga::PointF point{1.f, 2.f, 3.f, 1.F};
//
//    const auto a = point * point;
//    EXPECT_EQ(a[Scalar], -1);
//
//    tiny_pga::PlaneF plane{1.f, 2.f, 3.f, 4.f};
//    const auto b = plane * plane;
//    EXPECT_EQ(b[Scalar], 2 * 2 + 3 * 3 + 4 * 4);
//}
//
//TEST(BasicTest, RotorTest)
//{
//    tiny_pga::RotorF rotor{1.F, 0.F, 0.F, 0.F};
//    tiny_pga::PointF point{ 1.f, 2.f, 3.f, 1.F};
//
//    const auto a = rotor * point * ~rotor;
//
//    EXPECT_NEAR(a.scalar(), 0.F, 1e-5);
//    EXPECT_NEAR(a.e0(), 0.F, 1e-5);
//
//    EXPECT_NEAR(a.e021(), point.e021(), 1e-5);
//    EXPECT_NEAR(a.e013(), point.e013(), 1e-5);
//    EXPECT_NEAR(a.e032(), point.e032(), 1e-5);
//}
//
//TEST(BasicTest, ReverseTest)
//{
//    tiny_pga::RotorF rotor{11.F, 12.F, 13.F, 14.F};
//    tiny_pga::PointF point{1.f, 2.f, 3.f, 4.F};
//
//    tiny_pga::RotorF r = ~rotor;
//    tiny_pga::PointF p = ~point;
//
//    EXPECT_NEAR(r.e021(), -rotor.e021(), 1e-5);
//    EXPECT_NEAR(r.e013(), -rotor.e013(), 1e-5);
//    EXPECT_NEAR(r.e032(), -rotor.e032(), 1e-5);
//    EXPECT_NEAR(r.e123(), -rotor.e123(), 1e-5);
//
//    EXPECT_NEAR(p.e021(), -point.e021(), 1e-5);
//    EXPECT_NEAR(p.e013(), -point.e013(), 1e-5);
//    EXPECT_NEAR(p.e032(), -point.e032(), 1e-5);
//    EXPECT_NEAR(p.e123(), -point.e123(), 1e-5);
//}


//TEST(BasicTest, SandwichRotorTest)
//{
//    tiny_pga::RotorF rotor{{}, {1.F, 0.F, 0.F, 0.F}, {}};
//    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 1.F}};
//
//    tiny_pga::PointF a = point.sandwich(rotor);
//
//    EXPECT_NEAR(a.e021(), point.e021(), 1e-5);
//    EXPECT_NEAR(a.e013(), point.e013(), 1e-5);
//    EXPECT_NEAR(a.e032(), point.e032(), 1e-5);
//    EXPECT_NEAR(a.e123(), point.e123(), 1e-5);
//}


