//
// Created by Sergey Smirnov on 6.12.2020.
//

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

  constexpr Elems elem = static_cast<Elems>(elems::BitValues::kScalar);
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

TEST(BasicTest, MultiplicationAutoTest)
{
  tiny_pga::Plane plane{{1.f, 2.f, 3.f, 4.f}};

  tiny_pga::Point point{{}, {}, {}, {1.f, 2.f, 3.f, 0.F}};

  auto a = plane * point;

  tiny_pga::Rotor t {{}, {0.f, 1.f, 0.f, 0.f}};

  //tiny_pga::Plane other_plane = (t * plane) * ~t;

}