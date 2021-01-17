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

class PgaTest : public ::testing::Test
{
  public:
    PgaTest() = default;
};

TEST_F(PgaTest, BasicElemTest)
{
    constexpr Elems elem_1 = static_cast<Elems>(elems::BitValues::kScalar);
    Multivector<elem_1> mv{};

    EXPECT_TRUE(elems::has_scalar(static_cast<Elems>(elems::BitValues::kScalar)));
    EXPECT_TRUE(elems::has_scalar(Multivector<elem_1>::Elements));

    mv.scalar() = 15.f;
    EXPECT_EQ(mv.scalar(), 15.F);

    Multivector<elem_1> mv1{};
    Multivector<elem_1> mv2{};
    Elems new_elem = elems::multiplication(elem_1, elem_1);
    EXPECT_EQ(new_elem, elem_1);
    auto mv3 = mv1 * mv2;
    EXPECT_NO_THROW(mv3.scalar());
    mv3.scalar() = 1.F;
    EXPECT_EQ(mv3.scalar(), 1.F);
    mv3.scalar() = 2.F;
    EXPECT_EQ(mv3.scalar(), 2.F);

    constexpr Elems elem_e0 = static_cast<Elems>(elems::BitValues::kE0);

    Multivector<elem_e0> mv4{};
    Multivector<elem_e0> mv5 = mv1 * mv4;
    Multivector<elem_e0> mv6 = mv4 * mv1;
    Multivector<0> out0 = Multivector<elem_e0>{} * Multivector<elem_e0>{};

    constexpr Elems elem_e1 = static_cast<Elems>(elems::BitValues::kE1);
    Multivector<elem_1> out = Multivector<elem_e1>{} * Multivector<elem_e1>{};
}

TEST(BasicTest, SquaringTest)
{
    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 1.F}};

    const auto a = point * point;
    EXPECT_EQ(a.scalar(), -1);

    tiny_pga::PlaneF plane{{1.f, 2.f, 3.f, 4.f}};
    const auto b = plane * plane;
    EXPECT_EQ(b.scalar(), 2 * 2 + 3 * 3 + 4 * 4);
}

TEST(BasicTest, RotorTest)
{
    tiny_pga::RotorF rotor{{}, {1.F, 0.F, 0.F, 0.F}, {}};
    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 1.F}};

    const auto a = rotor * point * ~rotor;

    EXPECT_NEAR(a.scalar(), 0.F, 1e-5);
    EXPECT_NEAR(a.e0(), 0.F, 1e-5);

    EXPECT_NEAR(a.e021(), point.e021(), 1e-5);
    EXPECT_NEAR(a.e013(), point.e013(), 1e-5);
    EXPECT_NEAR(a.e032(), point.e032(), 1e-5);
}

TEST(BasicTest, ReverseTest)
{
    tiny_pga::RotorF rotor{{1.F, 0.F, 0.F, 0.F}, {11.F, 12.F, 13.F, 14.F}, {}};
    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 4.F}};

    tiny_pga::RotorF r = ~rotor;
    tiny_pga::PointF p = ~point;

    EXPECT_NEAR(r.e021(), -rotor.e021(), 1e-5);
    EXPECT_NEAR(r.e013(), -rotor.e013(), 1e-5);
    EXPECT_NEAR(r.e032(), -rotor.e032(), 1e-5);
    EXPECT_NEAR(r.e123(), -rotor.e123(), 1e-5);

    EXPECT_NEAR(p.e021(), -point.e021(), 1e-5);
    EXPECT_NEAR(p.e013(), -point.e013(), 1e-5);
    EXPECT_NEAR(p.e032(), -point.e032(), 1e-5);
    EXPECT_NEAR(p.e123(), -point.e123(), 1e-5);
}

//
TEST(BasicTest, SandwichRotorTest)
{
    tiny_pga::RotorF rotor{{}, {1.F, 0.F, 0.F, 0.F}, {}};
    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 1.F}};

    tiny_pga::PointF a = point.sandwich(rotor);

    EXPECT_NEAR(a.e021(), point.e021(), 1e-5);
    EXPECT_NEAR(a.e013(), point.e013(), 1e-5);
    EXPECT_NEAR(a.e032(), point.e032(), 1e-5);
    EXPECT_NEAR(a.e123(), point.e123(), 1e-5);
}