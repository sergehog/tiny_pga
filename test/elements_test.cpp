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

/// Converts tuple of 16 bools into Elems value
Elems elements(
    const std::tuple<bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool>
        elems)
{
    return elems::elements(std::get<0>(elems),
                           std::get<1>(elems),
                           std::get<2>(elems),
                           std::get<3>(elems),
                           std::get<4>(elems),
                           std::get<5>(elems),
                           std::get<6>(elems),
                           std::get<7>(elems),
                           std::get<8>(elems),
                           std::get<9>(elems),
                           std::get<10>(elems),
                           std::get<11>(elems),
                           std::get<12>(elems),
                           std::get<13>(elems),
                           std::get<14>(elems),
                           std::get<15>(elems));
}

/// Tests elems::geometric_product() function
TEST(BasicTest, GeometricProductElementsTest)
{
    const Elems ScalarElems = static_cast<Elems>(elems::Values::kScalar);
    const Elems OutScalarElems = elems::geometric_product(ScalarElems, ScalarElems);
    EXPECT_EQ(OutScalarElems, ScalarElems);

    const Elems ComplexElems = static_cast<Elems>(elems::Values::kScalar) | static_cast<Elems>(elems::Values::kE12);
    const Elems OutComplexElems = elems::geometric_product(ComplexElems, ComplexElems);
    EXPECT_EQ(OutComplexElems, ComplexElems);

    const Elems DualElems = static_cast<Elems>(elems::Values::kScalar) | static_cast<Elems>(elems::Values::kE0);
    const Elems OutDualElems = elems::geometric_product(DualElems, DualElems);
    EXPECT_EQ(OutDualElems, DualElems);

    // Multiplication of 2 rotors still a rotor
    const Elems OutRotorElems = elems::geometric_product(elems::RotorElems, elems::RotorElems);
    EXPECT_EQ(OutRotorElems, elems::RotorElems);

    // Multiplication of 2 translators still a translators
    const Elems OutTranslatorElems = elems::geometric_product(elems::TranslatorElems, elems::TranslatorElems);
    EXPECT_EQ(OutTranslatorElems, elems::TranslatorElems);

    // Multiplication of 2 motors still a motor
    const Elems OutMotorElems = elems::geometric_product(elems::MotorElems, elems::MotorElems);
    EXPECT_EQ(OutMotorElems, elems::MotorElems);
}

/// Test class for checking all possible element combinations
class PgaElementsTest
    : public ::testing::TestWithParam<
          std::tuple<bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool>>
{
  public:
    PgaElementsTest() = default;
};

/// Tests elems::elements() and family of "elems::has_ELEMENT()" functions, verifies that Elems
/// datatype is sufficient
TEST_P(PgaElementsTest, BasicElementsTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const Elems elem = elements(GetParam());

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

/// Certain combinations has to square to other certain elements
TEST_P(PgaElementsTest, GeometricProductSquaringTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const bool is_scalar_expected =
        test_scalar || test_e1 || test_e2 || test_e3 || test_e12 || test_e31 || test_e23 || test_e123;
    const bool is_e0_expected = (test_scalar && test_e0) || (test_e1 && test_e01) || (test_e2 && test_e02) ||
                                (test_e3 && test_e03) || (test_e12 && test_e021) || (test_e31 && test_e013) ||
                                (test_e23 && test_e032) || (test_e123 && test_e0123);

    const bool is_e1_expected =
        (test_scalar && test_e1) || (test_e2 && test_e12) || (test_e3 && test_e31) || (test_e23 && test_e123);
    const bool is_e2_expected =
        (test_scalar && test_e2) || (test_e1 && test_e12) || (test_e3 && test_e23) || (test_e31 && test_e123);
    const bool is_e3_expected =
        (test_scalar && test_e3) || (test_e1 && test_e31) || (test_e2 && test_e23) || (test_e12 && test_e123);

    const bool is_e01_expected = (test_scalar && test_e01) || (test_e0 && test_e1) || (test_e2 && test_e021) ||
                                 (test_e3 && test_e013) || (test_e02 && test_e12) || (test_e03 && test_e31) ||
                                 (test_e23 && test_e0123) || (test_e032 && test_e123);
    const bool is_e02_expected = (test_scalar && test_e02) || (test_e0 && test_e2) || (test_e1 && test_e021) ||
                                 (test_e3 && test_e032) || (test_e01 && test_e12) || (test_e03 && test_e23) ||
                                 (test_e31 && test_e0123) || (test_e013 && test_e123);
    const bool is_e03_expected = (test_scalar && test_e03) || (test_e0 && test_e3) || (test_e1 && test_e013) ||
                                 (test_e2 && test_e032) || (test_e01 && test_e31) || (test_e02 && test_e23) ||
                                 (test_e12 && test_e0123) || (test_e021 && test_e123);

    const bool is_e12_expected =
        (test_scalar && test_e12) || (test_e1 && test_e2) || (test_e3 && test_e123) || (test_e31 && test_e23);
    const bool is_e31_expected =
        (test_scalar && test_e31) || (test_e1 && test_e3) || (test_e2 && test_e123) || (test_e12 && test_e23);
    const bool is_e23_expected =
        (test_scalar && test_e23) || (test_e2 && test_e3) || (test_e1 && test_e123) || (test_e31 && test_e12);

    const bool is_e021_expected = (test_scalar && test_e021) || (test_e0 && test_e12) || (test_e1 && test_e02) ||
                                  (test_e2 && test_e01) || (test_e3 && test_e0123) || (test_e03 && test_e123) ||
                                  (test_e31 && test_e032) || (test_e23 && test_e013);
    const bool is_e013_expected = (test_scalar && test_e013) || (test_e0 && test_e31) || (test_e1 && test_e03) ||
                                  (test_e2 && test_e0123) || (test_e3 && test_e01) || (test_e02 && test_e123) ||
                                  (test_e12 && test_e032) || (test_e23 && test_e021);
    const bool is_e032_expected = (test_scalar && test_e032) || (test_e0 && test_e23) || (test_e1 && test_e0123) ||
                                  (test_e2 && test_e03) || (test_e3 && test_e02) || (test_e01 && test_e123) ||
                                  (test_e12 && test_e013) || (test_e31 && test_e021);

    const bool is_e123_expected =
        (test_scalar && test_e123) || (test_e1 && test_e23) || (test_e2 && test_e31) || (test_e3 && test_e12);

    const bool is_e0123_expected = (test_scalar && test_e0123) || (test_e0 && test_e123) || (test_e1 && test_e032) ||
                                   (test_e2 && test_e013) || (test_e3 && test_e021) || (test_e01 && test_e23) ||
                                   (test_e02 && test_e31) || (test_e03 && test_e12);

    const Elems elem = elems::elements(test_scalar,
                                       test_e0,
                                       test_e1,
                                       test_e2,
                                       test_e3,
                                       test_e01,
                                       test_e02,
                                       test_e03,
                                       test_e12,
                                       test_e31,
                                       test_e23,
                                       test_e021,
                                       test_e013,
                                       test_e032,
                                       test_e123,
                                       test_e0123);

    const Elems resulting_elems = elems::geometric_product(elem, elem);

    EXPECT_EQ(elems::has_scalar(resulting_elems), is_scalar_expected);
    EXPECT_EQ(elems::has_e0(resulting_elems), is_e0_expected);
    EXPECT_EQ(elems::has_e1(resulting_elems), is_e1_expected);
    EXPECT_EQ(elems::has_e2(resulting_elems), is_e2_expected);
    EXPECT_EQ(elems::has_e3(resulting_elems), is_e3_expected);
    EXPECT_EQ(elems::has_e01(resulting_elems), is_e01_expected);
    EXPECT_EQ(elems::has_e02(resulting_elems), is_e02_expected);
    EXPECT_EQ(elems::has_e03(resulting_elems), is_e03_expected);
    EXPECT_EQ(elems::has_e12(resulting_elems), is_e12_expected);
    EXPECT_EQ(elems::has_e31(resulting_elems), is_e31_expected);
    EXPECT_EQ(elems::has_e23(resulting_elems), is_e23_expected);
    EXPECT_EQ(elems::has_e021(resulting_elems), is_e021_expected);
    EXPECT_EQ(elems::has_e013(resulting_elems), is_e013_expected);
    EXPECT_EQ(elems::has_e032(resulting_elems), is_e032_expected);
    EXPECT_EQ(elems::has_e123(resulting_elems), is_e123_expected);
    EXPECT_EQ(elems::has_e0123(resulting_elems), is_e0123_expected);
}

/// Certain combinations has to square to other certain elements
TEST_P(PgaElementsTest, InnerProductSquaringTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const bool is_scalar_expected =
        test_scalar || test_e1 || test_e2 || test_e3 || test_e12 || test_e31 || test_e23 || test_e123;
    const bool is_e0_expected = (test_scalar && test_e0) || (test_e1 && test_e01) || (test_e2 && test_e02) ||
                                (test_e3 && test_e03) || (test_e12 && test_e021) || (test_e31 && test_e013) ||
                                (test_e23 && test_e032) || (test_e123 && test_e0123);
    const bool is_e1_expected =
        (test_scalar && test_e1) || (test_e2 && test_e12) || (test_e3 && test_e31) || (test_e23 && test_e123);
    const bool is_e2_expected =
        (test_scalar && test_e2) || (test_e1 && test_e12) || (test_e3 && test_e23) || (test_e31 && test_e123);
    const bool is_e3_expected =
        (test_scalar && test_e3) || (test_e1 && test_e31) || (test_e2 && test_e23) || (test_e12 && test_e123);
    const bool is_e01_expected =
        (test_scalar && test_e01) || (test_e2 && test_e021) || (test_e3 && test_e013) || (test_e23 && test_e0123);
    const bool is_e02_expected =
        (test_scalar && test_e02) || (test_e1 && test_e021) || (test_e3 && test_e032) || (test_e31 && test_e0123);
    const bool is_e03_expected =
        (test_scalar && test_e03) || (test_e1 && test_e013) || (test_e2 && test_e032) || (test_e12 && test_e0123);
    const bool is_e12_expected = (test_scalar && test_e12) || (test_e3 && test_e123);
    const bool is_e31_expected = (test_scalar && test_e31) || (test_e2 && test_e123);
    const bool is_e23_expected = (test_scalar && test_e23) || (test_e1 && test_e123);
    const bool is_e021_expected = (test_scalar && test_e021) || (test_e3 && test_e0123);
    const bool is_e013_expected = (test_scalar && test_e013) || (test_e2 && test_e0123);
    const bool is_e032_expected = (test_scalar && test_e032) || (test_e1 && test_e0123);
    const bool is_e123_expected = (test_scalar && test_e123);
    const bool is_e0123_expected = (test_scalar && test_e0123);

    const Elems elem = elems::elements(test_scalar,
                                       test_e0,
                                       test_e1,
                                       test_e2,
                                       test_e3,
                                       test_e01,
                                       test_e02,
                                       test_e03,
                                       test_e12,
                                       test_e31,
                                       test_e23,
                                       test_e021,
                                       test_e013,
                                       test_e032,
                                       test_e123,
                                       test_e0123);

    const Elems resulting_elems = elems::inner_product(elem, elem);
    EXPECT_EQ(elems::has_scalar(resulting_elems), is_scalar_expected);
    EXPECT_EQ(elems::has_e0(resulting_elems), is_e0_expected);
    EXPECT_EQ(elems::has_e1(resulting_elems), is_e1_expected);
    EXPECT_EQ(elems::has_e2(resulting_elems), is_e2_expected);
    EXPECT_EQ(elems::has_e3(resulting_elems), is_e3_expected);
    EXPECT_EQ(elems::has_e01(resulting_elems), is_e01_expected);
    EXPECT_EQ(elems::has_e02(resulting_elems), is_e02_expected);
    EXPECT_EQ(elems::has_e03(resulting_elems), is_e03_expected);
    EXPECT_EQ(elems::has_e12(resulting_elems), is_e12_expected);
    EXPECT_EQ(elems::has_e31(resulting_elems), is_e31_expected);
    EXPECT_EQ(elems::has_e23(resulting_elems), is_e23_expected);
    EXPECT_EQ(elems::has_e021(resulting_elems), is_e021_expected);
    EXPECT_EQ(elems::has_e013(resulting_elems), is_e013_expected);
    EXPECT_EQ(elems::has_e032(resulting_elems), is_e032_expected);
    EXPECT_EQ(elems::has_e123(resulting_elems), is_e123_expected);
    EXPECT_EQ(elems::has_e0123(resulting_elems), is_e0123_expected);
}

/// Certain combinations has to square to other certain elements
TEST_P(PgaElementsTest, OuterProductSquaringTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const bool is_scalar_expected = test_scalar;
    const bool is_e0_expected = (test_scalar && test_e0);
    const bool is_e1_expected = (test_scalar && test_e1);
    const bool is_e2_expected = (test_scalar && test_e2);
    const bool is_e3_expected = (test_scalar && test_e3);
    const bool is_e01_expected = (test_scalar && test_e01) || (test_e0 && test_e1);
    const bool is_e02_expected = (test_scalar && test_e02) || (test_e0 && test_e2);
    const bool is_e03_expected = (test_scalar && test_e03) || (test_e0 && test_e3);
    const bool is_e12_expected = (test_scalar && test_e12) || (test_e1 && test_e2);
    const bool is_e31_expected = (test_scalar && test_e31) || (test_e1 && test_e3);
    const bool is_e23_expected = (test_scalar && test_e23) || (test_e2 && test_e3);
    const bool is_e021_expected =
        (test_scalar && test_e021) || (test_e0 && test_e12) || (test_e1 && test_e02) || (test_e2 && test_e01);
    const bool is_e013_expected =
        (test_scalar && test_e013) || (test_e0 && test_e31) || (test_e1 && test_e03) || (test_e3 && test_e01);
    const bool is_e032_expected =
        (test_scalar && test_e032) || (test_e0 && test_e23) || (test_e2 && test_e03) || (test_e3 && test_e02);
    const bool is_e123_expected =
        (test_scalar && test_e123) || (test_e1 && test_e23) || (test_e2 && test_e31) || (test_e3 && test_e12);
    const bool is_e0123_expected = (test_scalar && test_e0123) || (test_e0 && test_e123) || (test_e1 && test_e032) ||
                                   (test_e2 && test_e013) || (test_e3 && test_e021) || (test_e01 && test_e23) ||
                                   (test_e02 && test_e31) || (test_e03 && test_e12);
    const Elems elem = elems::elements(test_scalar,
                                       test_e0,
                                       test_e1,
                                       test_e2,
                                       test_e3,
                                       test_e01,
                                       test_e02,
                                       test_e03,
                                       test_e12,
                                       test_e31,
                                       test_e23,
                                       test_e021,
                                       test_e013,
                                       test_e032,
                                       test_e123,
                                       test_e0123);

    const Elems resulting_elems = elems::outer_product(elem, elem);
    EXPECT_EQ(elems::has_scalar(resulting_elems), is_scalar_expected);
    EXPECT_EQ(elems::has_e0(resulting_elems), is_e0_expected);
    EXPECT_EQ(elems::has_e1(resulting_elems), is_e1_expected);
    EXPECT_EQ(elems::has_e2(resulting_elems), is_e2_expected);
    EXPECT_EQ(elems::has_e3(resulting_elems), is_e3_expected);
    EXPECT_EQ(elems::has_e01(resulting_elems), is_e01_expected);
    EXPECT_EQ(elems::has_e02(resulting_elems), is_e02_expected);
    EXPECT_EQ(elems::has_e03(resulting_elems), is_e03_expected);
    EXPECT_EQ(elems::has_e12(resulting_elems), is_e12_expected);
    EXPECT_EQ(elems::has_e31(resulting_elems), is_e31_expected);
    EXPECT_EQ(elems::has_e23(resulting_elems), is_e23_expected);
    EXPECT_EQ(elems::has_e021(resulting_elems), is_e021_expected);
    EXPECT_EQ(elems::has_e013(resulting_elems), is_e013_expected);
    EXPECT_EQ(elems::has_e032(resulting_elems), is_e032_expected);
    EXPECT_EQ(elems::has_e123(resulting_elems), is_e123_expected);
    EXPECT_EQ(elems::has_e0123(resulting_elems), is_e0123_expected);
}

// TEST_P(PgaElementsTest, GeometricProductIsSumOfInnerAndOuterTest)
//{
//    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
//        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
//    std::tie(test_scalar,
//             test_e0,
//             test_e1,
//             test_e2,
//             test_e3,
//             test_e01,
//             test_e02,
//             test_e03,
//             test_e12,
//             test_e31,
//             test_e23,
//             test_e021,
//             test_e013,
//             test_e032,
//             test_e123,
//             test_e0123) = GetParam();
//
//    const Elems elem = elems::elements(test_scalar,
//                                       test_e0,
//                                       test_e1,
//                                       test_e2,
//                                       test_e3,
//                                       test_e01,
//                                       test_e02,
//                                       test_e03,
//                                       test_e12,
//                                       test_e31,
//                                       test_e23,
//                                       test_e021,
//                                       test_e013,
//                                       test_e032,
//                                       test_e123,
//                                       test_e0123);
//
//    const Elems gp_elems = elems::geometric_product(elem, elem);
//    const Elems ip_elems = elems::inner_product(elem, elem);
//    const Elems op_elems = elems::outer_product(elem, elem);
//    const Elems test_elems = ip_elems | op_elems;
//
//    EXPECT_EQ(gp_elems, test_elems);
//}

/// Dualization of elements two times must result in the same elements as initial
TEST_P(PgaElementsTest, DualElementsTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const Elems elems = elems::elements(test_scalar,
                                        test_e0,
                                        test_e1,
                                        test_e2,
                                        test_e3,
                                        test_e01,
                                        test_e02,
                                        test_e03,
                                        test_e12,
                                        test_e31,
                                        test_e23,
                                        test_e021,
                                        test_e013,
                                        test_e032,
                                        test_e123,
                                        test_e0123);

    const Elems dual_elems = elems::dual(elems);

    // some "symmetric" element combinations might be dual of each other
    // here we check that at least some duals are not the same as original
    if (elems != 0U && elems::count(elems) % 2)
    {
        EXPECT_NE(elems, dual_elems);
    }

    // however, double dual must be the same, as original
    const Elems double_dual_elems = elems::dual(dual_elems);
    EXPECT_EQ(elems, double_dual_elems);
}

TEST_P(PgaElementsTest, ElemsAdditionTest)
{
    bool test_scalar{}, test_e0{}, test_e1{}, test_e2{}, test_e3{}, test_e01{}, test_e02{}, test_e03{}, test_e12{},
        test_e31{}, test_e23{}, test_e021{}, test_e013{}, test_e032{}, test_e123{}, test_e0123{};
    std::tie(test_scalar,
             test_e0,
             test_e1,
             test_e2,
             test_e3,
             test_e01,
             test_e02,
             test_e03,
             test_e12,
             test_e31,
             test_e23,
             test_e021,
             test_e013,
             test_e032,
             test_e123,
             test_e0123) = GetParam();

    const Elems elems = elems::elements(test_scalar,
                                        test_e0,
                                        test_e1,
                                        test_e2,
                                        test_e3,
                                        test_e01,
                                        test_e02,
                                        test_e03,
                                        test_e12,
                                        test_e31,
                                        test_e23,
                                        test_e021,
                                        test_e013,
                                        test_e032,
                                        test_e123,
                                        test_e0123);

    const Elems empty = 0U;

    const Elems elems1 = elems::addition(elems, empty);
    const Elems elems2 = elems::addition(empty, elems);

    const Elems scalar = static_cast<Elems>(elems::Values::kScalar);
    const Elems elems1s = elems::addition(elems, scalar);
    const Elems elems2s = elems::addition(scalar, elems);

    EXPECT_EQ(elems, elems1);
    EXPECT_EQ(elems, elems2);
    EXPECT_EQ(elems1s, elems2s);
    EXPECT_TRUE(elems::has_scalar(elems1s));
}

INSTANTIATE_TEST_CASE_P(AllElements,
                        PgaElementsTest,
                        testing::Combine(testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool(),
                                         testing::Bool()));
