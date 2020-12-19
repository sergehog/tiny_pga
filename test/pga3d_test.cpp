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

    const auto M1 = motor_from_3_points_pairs({A, B, C}, {A1, B1, C1});
    const auto M2 = motor_from_3_points_pairs({A, C, B}, {A1, C1, B1});

    EXPECT_NEAR(M1[0], M2[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M2[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M2[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M2[kE03], 1E-7);
    EXPECT_NEAR(M1[kE12], M2[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M2[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M2[kE31], 1E-7);

    const auto M3 = motor_from_3_points_pairs({C, A, B}, {C1, A1, B1});

    EXPECT_NEAR(M1[0], M3[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M3[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M3[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M3[kE03], 2E-7);
    EXPECT_NEAR(M1[kE12], M3[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M3[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M3[kE31], 1E-7);

    const auto M4 = motor_from_3_points_pairs({B, A, C}, {B1, A1, C1});

    EXPECT_NEAR(M1[0], M4[0], 1E-7);
    EXPECT_NEAR(M1[kE01], M4[kE01], 1E-7);
    EXPECT_NEAR(M1[kE02], M4[kE02], 1E-7);
    EXPECT_NEAR(M1[kE03], M4[kE03], 2E-7);
    EXPECT_NEAR(M1[kE12], M4[kE12], 1E-7);
    EXPECT_NEAR(M1[kE23], M4[kE23], 1E-7);
    EXPECT_NEAR(M1[kE31], M4[kE31], 1E-7);
}
