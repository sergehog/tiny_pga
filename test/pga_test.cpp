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

// TEST(BasicTest, SandwichAutoTest)
//{
//    tiny_pga::PointF point{{}, {}, {}, {1.f, 2.f, 3.f, 0.F}};
//
//    tiny_pga::TranslatorF t {{},{5.f, 7.f, 9.f, 11.F}};
//
//    tiny_pga::PointF p2 = point.sandwich(t);
//    EXPECT_EQ(p2.e01(), 11.F);
//}

/*
template <Elems elems>
struct AllPgaTest
{
    static void RunTest()
    {
        tiny_pga::Multivector<elems> pga {};
        if(elems::has_scalar(elems))
        {
            pga.scalar() = 1.F;
        }
        if(elems::has_e0(elems))
        {
            pga.e0() = 1.F;
        }
        if(elems::has_e1(elems))
        {
            pga.e1() = 1.F;
        }
        if(elems::has_e2(elems))
        {
            pga.e2() = 1.F;
        }
        if(elems::has_e3(elems))
        {
            pga.e3() = 1.F;
        }
        if(elems::has_e01(elems))
        {
            pga.e01() = 1.F;
        }
        if(elems::has_e02(elems))
        {
            pga.e02() = 1.F;
        }
        if(elems::has_e03(elems))
        {
            pga.e03() = 1.F;
        }
        if(elems::has_e31(elems))
        {
            pga.e31() = 1.F;
        }
        if(elems::has_e23(elems))
        {
            pga.e23() = 1.F;
        }
        if(elems::has_e12(elems))
        {
            pga.e12() = 1.F;
        }
        if(elems::has_e021(elems))
        {
            pga.e021() = 1.F;
        }
        if(elems::has_e013(elems))
        {
            pga.e013() = 1.F;
        }
        if(elems::has_e032(elems))
        {
            pga.e032() = 1.F;
        }
        if(elems::has_e123(elems))
        {
            pga.e123() = 1.F;
        }
        if(elems::has_e0123(elems))
        {
            pga.e0123() = 1.F;
        }

        constexpr Elems res_elems = elems::multiplication(elems, elems);
        tiny_pga::Multivector<res_elems> result = pga * pga;

        if(elems::has_scalar(res_elems))
        {
            float scalar = 0.F;
            scalar += elems::has_scalar(elems);
            scalar += elems::has_e1(elems);
            scalar += elems::has_e2(elems);
            scalar += elems::has_e3(elems);
            scalar -= elems::has_e12(elems);
            scalar -= elems::has_e31(elems);
            scalar -= elems::has_e23(elems);
            scalar -= elems::has_e123(elems);

            EXPECT_EQ(pga.scalar(), scalar);
        }



        if(elems < std::numeric_limits<Elems>::max())
        {
            AllPgaTest<elems+1>::RunTest();
        }
    }

};
TEST(BasicTest, MultivectorSquaringTest)
{
    AllPgaTest<0>::RunTest();
}
*/