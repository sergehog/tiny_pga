// This code was originally downloaded from bivector.net
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
#include "../pga3d.h"
#include <gtest/gtest.h>

using namespace tiny_pga;
using namespace float_basis;

TEST(PGA3DTest, BasicTest)
{
    // Elements of the even subalgebra (scalar + bivector + pss) of unit length are motors
    PGA3D<> rot = rotor(PI / 2.0f, e1 * e2);

    // The outer product ^ is the MEET. Here we intersect the yz (x=0) and xz (y=0) planes.
    PGA3D<> ax_z = e1 ^ e2;

    // line and plane meet in point. We intersect the line along the z-axis (x=0,y=0) with the xy (z=0) plane.
    PGA3D<> orig = ax_z ^ e3;

    // We can also easily create points and join them into a line using the regressive (vee, &) product.
    PGA3D<> px = point(1.f, 0.f, 0.f);
    PGA3D<> line = orig & px;

    // Lets also create the plane with equation 2x + z - 3 = 0
    PGA3D<> p = plane(2.f, 0.f, 1.f, -3.f);

    // rotations work on all elements
    PGA3D<> rotated_plane = rot * p * ~rot;
    PGA3D<> rotated_line = rot * line * ~rot;
    PGA3D<> rotated_point = rot * px * ~rot;

    // See the 3D PGA Cheat sheet for a huge collection of useful formulas
    PGA3D<> point_on_plane = (p | px) * p;

    // Some output.
    printf("a point       : ");
    log(px);
    printf("a line        : ");
    log(line);
    printf("a plane       : ");
    log(p);
    printf("a rotor       : ");
    log(rot);
    printf("rotated line  : ");
    log(rotated_line);
    printf("rotated point : ");
    log(rotated_point);
    printf("rotated plane : ");
    log(rotated_plane);
    printf("point on plane: ");
    log(point_on_plane.normalized());
    printf("point on torus: ");
    log(point_on_torus(0.0f, 0.0f));
    log(e0 - 1.0f);
    log(1.0f - e0);
}

TEST(PGA3DTest, ConstMotorEstimatorTest)
{
    PGA3D<> A = point(0.F, 0.F, 0.F);
    PGA3D<> B = point(1.F, 0.F, 0.F);
    PGA3D<> C = point(0.F, 1.F, 0.F);

    PGA3D<> A1 = point(1.F, 1.F, 1.F);
    PGA3D<> B1 = point(1.F, 2.F, 1.F);
    PGA3D<> C1 = point(1.F, 1.F, 2.F);

    // Va = sqrt(Ai/A);
    const auto Va = (1.F + (A1 * ~A)).normalized();
    log(Va, "Va");

    const auto Ba = Va * B * ~Va;
    log(Ba, "Ba");

    const auto Vb_squared = (A1 & B1) * ~(A1 & Ba);
    log(Vb_squared, "Vb_squared");

    const auto Vb = (1 + Vb_squared).normalized();  // Vb = sqrt(Vb_squared);
    log(Vb, "Vb");

    const auto Cba = Vb * Va * C * ~Va * ~Vb;
    log(Cba, "Cba");

    const auto Vc_squared = (A1 & B1 & C1) * ~(A1 & B1 & Cba);
    log(Vc_squared, "Vc_squared");

    const auto Vc = (1 + Vc_squared).normalized();  // Vc = sqrt(Vc_squared);
    log(Vc, "Vc");

    const auto V = Vc * Vb * Va;
    log(V, "V");

    const auto Ai = V * A * ~V - A1;
    log(Ai, "V*A*~V - A1");

    const auto Bi = V * B * ~V - B1;
    log(Bi, "V*B*~V - B1");

    const auto Ci = V * C * ~V - C1;
    log(Ci, "V*C*~V - C1");

    EXPECT_NEAR((Ai * Ai)[0], 0.F, 0.0001);
    EXPECT_NEAR((Bi * Bi)[0], 0.F, 0.0001);
    EXPECT_NEAR((Bi * Bi)[0], 0.F, 0.0001);
}

/// Estimated Motor shall be the same, despite order of given points
TEST(PGA3DTest, MotorEstimatorStabilityTest)
{
    PGA3D<> A = point(1.F, 1.F, 1.F);
    PGA3D<> B = point(2.F, 1.F, 1.F);
    PGA3D<> C = point(1.F, 2.F, 1.F);

    PGA3D<> A1 = point(-1.F, -1.F, -1.F);
    PGA3D<> B1 = point(-1.F, -2.F, -1.F);
    PGA3D<> C1 = point(-1.F, -1.F, -2.F);

    const auto M1 = motor_from_point_pairs({A, B, C}, {A1, B1, C1});
    const auto M2 = motor_from_point_pairs({A, C, B}, {A1, C1, B1});

    EXPECT_NEAR(M1[0], M2[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M2[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M2[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M2[kE03], 1E-7);
    EXPECT_NEAR(M1[kE12], M2[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M2[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M2[kE31], 1E-7);

    const auto M3 = motor_from_point_pairs({C, A, B}, {C1, A1, B1});

    EXPECT_NEAR(M1[0], M3[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M3[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M3[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M3[kE03], 2E-7);
    EXPECT_NEAR(M1[kE12], M3[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M3[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M3[kE31], 1E-7);

    const auto M4 = motor_from_point_pairs({B, A, C}, {B1, A1, C1});

    EXPECT_NEAR(M1[0], M4[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M4[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M4[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M4[kE03], 2E-7);
    EXPECT_NEAR(M1[kE12], M4[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M4[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M4[kE31], 1E-7);
}

TEST(PGA3DTest, How_AverageMotor_ComparesTo_MotorOfAverage_Test)
{
    const PGA3D<> A1 = point(1.F, 1.F, 1.F);
    const PGA3D<> A2 = point(2.F, 1.F, 1.F);
    const PGA3D<> A3 = point(2.F, 2.F, 1.F);

    const PGA3D<> B1 = point(-1.F, -1.F, -1.F);
    const PGA3D<> B2 = point(-1.F, -2.F, -1.F);
    const PGA3D<> B3 = point(-1.F, -2.F, -2.F);

    const PGA3D<> Aavg = (A1 + A2 + A3) * (1.F / 3.F);
    const PGA3D<> Bavg = (B1 + B2 + B3) * (1.f / 3.F);

    const PGA3D<> M1 = (B1 * ~A1).sqrt();
    const PGA3D<> M2 = (B2 * ~A2).sqrt();
    const PGA3D<> M3 = (B3 * ~A3).sqrt();

    // average motor
    const PGA3D<> M123avg = (M1 + M2 + M3) * (1.f / 3.F);

    // motor of average
    const PGA3D<> Mavg = (Bavg * ~Aavg).sqrt();

    log(M123avg, "M123avg");
    log(Mavg, "Mavg");
    EXPECT_NEAR(M123avg[0], Mavg[0], 1E-7);
    EXPECT_NEAR(M123avg[kE01], Mavg[kE01], 1E-7);
    EXPECT_NEAR(M123avg[kE02], Mavg[kE02], 1E-7);
    EXPECT_NEAR(M123avg[kE03], Mavg[kE03], 2E-7);
    EXPECT_NEAR(M123avg[kE12], Mavg[kE12], 1E-7);
    EXPECT_NEAR(M123avg[kE23], Mavg[kE23], 1E-7);
    EXPECT_NEAR(M123avg[kE31], Mavg[kE31], 1E-7);
}

template <typename ScalarType = float>
PGA3D<ScalarType> motor_from_3point_pairs(const std::array<PGA3D<ScalarType>, 3>& reference_points,
                                          const std::array<PGA3D<ScalarType>, 3>& target_points)
{
    const PGA3D<>& A1 = reference_points[0];
    const PGA3D<>& A2 = reference_points[1];
    const PGA3D<>& A3 = reference_points[2];
    const PGA3D<> B1 = target_points[0];
    const PGA3D<> B2 = target_points[1];
    const PGA3D<> B3 = target_points[2];

    const PGA3D<> Aavg = (A1 + A2 + A3) * (1.F / 3.F);
    const PGA3D<> Bavg = (B1 + B2 + B3) * (1.f / 3.F);

    // motor of average for points
    const PGA3D<> Mavg = (Bavg * ~Aavg).sqrt();

    // "translated" points
    const PGA3D<> A1a = Mavg * A1 * ~Mavg;
    const PGA3D<> A2a = Mavg * A2 * ~Mavg;
    const PGA3D<> A3a = Mavg * A3 * ~Mavg;

    // motors of lines
    const PGA3D<> M1a = ((Bavg & B1) * ~(Bavg & A1a)).sqrt();
    const PGA3D<> M2a = ((Bavg & B2) * ~(Bavg & A2a)).sqrt();
    const PGA3D<> M3a = ((Bavg & B3) * ~(Bavg & A3a)).sqrt();
    // average motor of lines
    const PGA3D<> M123a = (M1a + M2a + M3a) * (1.f / 3.F);

    // average line
    const auto Bline = ((Bavg & B1) + (Bavg & B2) + (Bavg & B3)) * (1.F / 3.F);

    // "translated" points again
    const PGA3D<> A1b = M123a * A1a * ~M123a;
    const PGA3D<> A2b = M123a * A2a * ~M123a;
    const PGA3D<> A3b = M123a * A3a * ~M123a;

    const PGA3D<> M1b = ((Bline & Bavg & B1) * ~(Bline & Bavg & A1b)).sqrt();
    const PGA3D<> M2b = ((Bline & Bavg & B2) * ~(Bline & Bavg & A2b)).sqrt();
    const PGA3D<> M3b = ((Bline & Bavg & B3) * ~(Bline & Bavg & A3b)).sqrt();
    // average motor of planes
    const PGA3D<> M123ab = (M1a + M2a + M3a) * (1.f / 3.F);

    return M123ab * M123a * Mavg;
}

template <typename ScalarType = float>
PGA3D<ScalarType> motor_from_many_point_pairs(const std::vector<PGA3D<ScalarType>>& reference_points,
                                              const std::vector<PGA3D<ScalarType>>& target_points)
{
    PGA3D<ScalarType> M1;
    std::size_t size = std::min(reference_points.size(), target_points.size());

    if (size == 0U)
    {
        return M1;
    }

    // Find centroids for reference and for target sets
    for (std::size_t i = 0U; i < size; i++)
    {
        M1 = M1 + target_points[i] * ~reference_points[i];
    }

    return M1.sqrt();
}

// template <typename ScalarType = float>
// PGA3D<ScalarType> motor_from_element_pairs(const std::vector<PGA3D<ScalarType>>& reference_elems,
//                                           const std::vector<PGA3D<ScalarType>>& target_elems)
//{
//    PGA3D<ScalarType> M1;
//    std::size_t size = std::min(reference_elems.size(), target_elems.size());
//
//    if (size == 0U)
//    {
//        return M1;
//    }
//
//    // Find first "average" transform (it behaves almost like finding centroid/average translation)
//    for (std::size_t i = 0U; i < size; i++)
//    {
//        // (A1 * ~A).sqrt();
//        auto Mi = (target_elems[i] * ~reference_elems[i]).sqrt();
//        M1 = M1 + Mi;
//    }
//    M1 = M1 * float(1.F / float(size));
//
//    // transfrom all reference elements to a new transform
//    std::vector<PGA3D<ScalarType>> reference_shifted;
//    for (std::size_t i = 0U; i < size; i++)
//    {
//
//        reference_shifted.push_back(M1 * reference_elems[i] * ~M1);
//    }
//
//    PGA3D<ScalarType> M2;
//    for (std::size_t i = 0U; i < size; i++)
//    {
//        //const auto Vb = ((A1 & B1) * ~(A1 & Ba)).sqrt();
//        auto Mi = ((M1 & target_elems[i]) * ~(M1 & reference_shifted[i])).sqrt();
//        M2 = M2 + Mi;
//    }
//    M2 = M2 * float(1.F / float(size));
//
////    // transfrom all reference elements to a new transform
////    for (std::size_t i = 0U; i < size; i++)
////    {
////        reference_shifted[i] = (M2 * reference_shifted[i] * ~M2);
////    }
////
////    PGA3D<ScalarType> M3;
////    for (std::size_t i = 0U; i < size; i++)
////    {
////        auto Mi = (target_elems[i] * ~reference_shifted[i]).sqrt();
////        M3 = M3 + Mi;
////    }
////    M3 = M3 * float(1.F / float(size));
//
//
//    return /*M3 **/ M2 * M1;
//}

// TEST(PGA3DTest, MotorEstimator3StabilityTest)
//{
//    PGA3D<> A = point(1.F, 1.F, 1.F);
//    PGA3D<> B = point(2.F, 1.F, 1.F);
//    PGA3D<> C = point(1.F, 2.F, 1.F);
//
//    PGA3D<> A1 = point(-1.F, -1.F, -1.F);
//    PGA3D<> B1 = point(-1.F, -2.F, -1.F);
//    PGA3D<> C1 = point(-1.F, -1.F, -2.F);
//
//    const auto M1 = motor_from_point_pairs({A, B, C}, {A1, B1, C1});
//    const auto M2 = motor_from_3point_pairs({A, C, B}, {A1, C1, B1});
//    log(M1, "M1");
//    log(M2, "M2");
//
//    EXPECT_NEAR(M1[0], M2[0], 1E-7);
//    EXPECT_NEAR(M1[kE01], M2[kE01], 1E-7);
//    EXPECT_NEAR(M1[kE02], M2[kE02], 1E-7);
//    EXPECT_NEAR(M1[kE03], M2[kE03], 1E-7);
//    EXPECT_NEAR(M1[kE12], M2[kE12], 1E-7);
//    EXPECT_NEAR(M1[kE23], M2[kE23], 1E-7);
//    EXPECT_NEAR(M1[kE31], M2[kE31], 1E-7);
//}
//
//
//
//
///// Estimated Motor shall be the same, despite order of given points
// TEST(PGA3DTest, MultipointMotorEstimatorTest)
//{
//    PGA3D<> A = point(1.F, 1.F, 1.F);
//    PGA3D<> B = point(2.F, 1.F, 1.F);
//    PGA3D<> C = point(2.F, 2.F, 1.F);
//    PGA3D<> D = point(1.F, 2.F, 1.F);
//
//
//    PGA3D<> A1 = point(-1.F, -1.F, -1.F);
//    PGA3D<> B1 = point(-1.F, -2.F, -1.F);
//    PGA3D<> C1 = point(-1.F, -2.F, -2.F);
//    PGA3D<> D1 = point(-1.F, -1.F, -2.F);
//
//    const auto M1 = motor_from_point_pairs({A, B, C}, {A1, B1, C1});
//    const auto M2 = motor_from_many_point_pairs({A, C, B}, {A1, C1, B1});
//
//    log(M1, "M1");
//    log(M2, "M2");
//
//    EXPECT_NEAR(M1[0], M2[0], 1E-7);
////    EXPECT_NEAR(M1[kE01], M2[kE01], 1E-7);
////    EXPECT_NEAR(M1[kE02], M2[kE02], 1E-7);
////    EXPECT_NEAR(M1[kE03], M2[kE03], 2E-7);
////    EXPECT_NEAR(M1[kE12], M2[kE12], 1E-7);
////    EXPECT_NEAR(M1[kE23], M2[kE23], 1E-7);
////    EXPECT_NEAR(M1[kE31], M2[kE31], 1E-7);
//
//}

// TEST(PGA3DTest, How_AverageMotor_ComparesTo_MotorOfAverage_ForLines_Test)
//{
//    const PGA3D<> A1 = point(1.F, 1.F, 1.F);
//    const PGA3D<> A2 = point(2.F, 1.F, 1.F);
//    const PGA3D<> A3 = point(2.F, 2.F, 1.F);
//
//    const PGA3D<> B1 = point(-1.F, -1.F, -1.F);
//    const PGA3D<> B2 = point(-1.F, -2.F, -1.F);
//    const PGA3D<> B3 = point(-1.F, -2.F, -2.F);
//
//    const PGA3D<> Aavg = (A1 + A2 + A3) * (1.F/3.F);
//    const PGA3D<> Bavg = (B1 + B2 + B3) * (1.f/3.F);
//
//    // motor of average for points
//    const PGA3D<> Mavg = (Bavg * ~Aavg).sqrt();
//
//    // "translated" points
//    const PGA3D<> A1a = Mavg * A1 * ~Mavg;
//    const PGA3D<> A2a = Mavg * A2 * ~Mavg;
//    const PGA3D<> A3a = Mavg * A3 * ~Mavg;
//
//    // motors of lines
//    const PGA3D<> M1a = ((Bavg & B1) * ~(Bavg & A1a)).sqrt();
//    const PGA3D<> M2a = ((Bavg & B2) * ~(Bavg & A2a)).sqrt();
//    const PGA3D<> M3a = ((Bavg & B3) * ~(Bavg & A3a)).sqrt();
//    // average motor of lines
//    const PGA3D<> M123a = (M1a + M2a + M3a)* (1.f/3.F);
//
//    // average lines
//    //const auto Aline = ((Bavg & A1a) + (Bavg & A2a) + (Bavg & A3a)) * (1.F/3.F);
//    const auto Bline = ((Bavg & B1) + (Bavg & B2) + (Bavg & B3)) * (1.F/3.F);
//
////    // motor of average lines
////    const PGA3D<> M_ABline = (Bline * ~Aline).sqrt();
////    log(M123a, "M123a");
////    log(M_ABline, "M_ABline");
////    EXPECT_NEAR(M123a[0], M_ABline[0], 1E-7);
//
//    // "translated" points again
//    const PGA3D<> A1b = M123a * A1a * ~M123a;
//    const PGA3D<> A2b = M123a * A2a * ~M123a;
//    const PGA3D<> A3b = M123a * A3a * ~M123a;
//
//    const PGA3D<> M1b = ((Bline & Bavg & B1) * ~(Bline & Bavg & A1b)).sqrt();
//    const PGA3D<> M2b = ((Bline & Bavg & B2) * ~(Bline & Bavg & A2b)).sqrt();
//    const PGA3D<> M3b = ((Bline & Bavg & B3) * ~(Bline & Bavg & A3b)).sqrt();
//    // average motor of planes
//    const PGA3D<> M123ab = (M1a + M2a + M3a)* (1.f/3.F);
//
//
//}