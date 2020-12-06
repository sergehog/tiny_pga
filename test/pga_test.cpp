//
// Created by Sergey Smirnov on 6.12.2020.
//

#include <gtest/gtest.h>
#include "../tiny_pga.h"


using namespace tiny_pga;

TEST(BasicTest, BasicElemTest)
{
  constexpr Elems elem = static_cast<Elems>(elems::BitValues::kScalar);
  Multivector<elem> mv{};

  EXPECT_TRUE(elems::has_scalar(static_cast<Elems>(elems::BitValues::kScalar)));
  EXPECT_TRUE(elems::has_scalar(Multivector<elem>::Elements));

  mv.scalar() = 15.f;
  EXPECT_EQ(mv.scalar(), 15.F);
}


TEST(BasicTest, MultiplicationElementsTest)
{
  constexpr Elems elem = static_cast<Elems>(elems::BitValues::kScalar);
  Multivector<elem> mv1 {};
  Multivector<elem> mv2 {};
  Elems new_elem = elems::multiplication(elem, elem);
  EXPECT_EQ(new_elem, elem);
  Multivector<elems::multiplication(elem, elem)> mv3;
  mv3 = mv1 * mv2;

  //auto mv3 = mv1 * mv2;
  //EXPECT_NO_THROW(mv3.scalar());
}
/*

TEST(BasicTest, MultiplicationTest)
{
  tiny_pga::Plane plane{1.f, 2.f, 3.f, 4.f};

  tiny_pga::Point point{{}, {}, {}, 1.f, 2.f, 3.f, 0.F};



}
*/
