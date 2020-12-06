//
// Created by Sergey Smirnov on 6.12.2020.
//

#include "../tiny_pga.h"
#include <gtest/gtest.h>

using namespace tiny_pga;

TEST(BasicTest, MultiplicationElementsTest)
{
  // after sandwich with a rotor plane must remain a plane
  Elems OutElems = elems::multiplication(
      elems::multiplication(elems::RotorElems, elems::PlaneElems), elems::RotorElems);
  EXPECT_EQ(OutElems, elems::PlaneElems);
}

/// Tests has_ELEMENT() family of functions
class PgaElementsTest
    : public ::testing::TestWithParam<std::tuple<bool, bool, bool, bool, bool, bool, bool, bool,
                                                 bool, bool, bool, bool, bool, bool, bool, bool>>
{
public:
  PgaElementsTest() = default;
};

TEST_P(PgaElementsTest, SimpleElementsTest)
{
  bool test_scalar, test_e0, test_e1, test_e2, test_e3, test_e01, test_e02, test_e03, test_e12,
      test_e31, test_e23, test_e021, test_e013, test_e032, test_e123, test_e0123;
  std::tie(test_scalar, test_e0, test_e1, test_e2, test_e3, test_e01, test_e02, test_e03, test_e12,
           test_e31, test_e23, test_e021, test_e013, test_e032, test_e123, test_e0123) = GetParam();

  const Elems elem = (test_scalar ? static_cast<Elems>(elems::BitValues::kScalar) : 0) |
                     (test_e0 ? static_cast<Elems>(elems::BitValues::kE0) : 0) |
                     (test_e1 ? static_cast<Elems>(elems::BitValues::kE1) : 0) |
                     (test_e2 ? static_cast<Elems>(elems::BitValues::kE2) : 0) |
                     (test_e3 ? static_cast<Elems>(elems::BitValues::kE3) : 0) |
                     (test_e01 ? static_cast<Elems>(elems::BitValues::kE01) : 0) |
                     (test_e02 ? static_cast<Elems>(elems::BitValues::kE02) : 0) |
                     (test_e03 ? static_cast<Elems>(elems::BitValues::kE03) : 0) |
                     (test_e12 ? static_cast<Elems>(elems::BitValues::kE12) : 0) |
                     (test_e31 ? static_cast<Elems>(elems::BitValues::kE31) : 0) |
                     (test_e23 ? static_cast<Elems>(elems::BitValues::kE23) : 0) |
                     (test_e021 ? static_cast<Elems>(elems::BitValues::kE021) : 0) |
                     (test_e013 ? static_cast<Elems>(elems::BitValues::kE013) : 0) |
                     (test_e032 ? static_cast<Elems>(elems::BitValues::kE032) : 0) |
                     (test_e123 ? static_cast<Elems>(elems::BitValues::kE123) : 0) |
                     (test_e0123 ? static_cast<Elems>(elems::BitValues::kE0123) : 0);

  EXPECT_EQ(elems::has_scalar(elem), test_scalar);
  EXPECT_EQ(elems::has_e0(elem), test_e0);
  EXPECT_EQ(elems::has_e1(elem), test_e1);
  EXPECT_EQ(elems::has_e2(elem), test_e2);
  EXPECT_EQ(elems::has_e3(elem), test_e3);
  EXPECT_EQ(elems::has_e01(elem), test_e01);
  EXPECT_EQ(elems::has_e02(elem), test_e02);
  EXPECT_EQ(elems::has_e03(elem), test_e03);
  EXPECT_EQ(elems::has_e12(elem), test_e12);
  EXPECT_EQ(elems::has_e31(elem), test_e31);
  EXPECT_EQ(elems::has_e23(elem), test_e23);
  EXPECT_EQ(elems::has_e021(elem), test_e021);
  EXPECT_EQ(elems::has_e013(elem), test_e013);
  EXPECT_EQ(elems::has_e032(elem), test_e032);
  EXPECT_EQ(elems::has_e123(elem), test_e123);
  EXPECT_EQ(elems::has_e0123(elem), test_e0123);
}

// INSTANTIATE_TEST_CASE_P(InstantiationName,
//                        PgaElementsTest,
//                        testing::Combine(
//                            testing::Bool(),testing::Bool(),testing::Bool(),testing::Bool(),
//                            testing::Bool(),testing::Bool(),testing::Bool(),testing::Bool(),
//                            testing::Bool(),testing::Bool(),testing::Bool(),testing::Bool(),
//                            testing::Bool(),testing::Bool(),testing::Bool(),testing::Bool())
//);
//
