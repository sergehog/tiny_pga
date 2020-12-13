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

TEST(PGA3DTest, BasicTest)
{
    const PGA3D<> e0(1.F, kE0);
    const PGA3D<> e1(1.F, kE1);
    const PGA3D<> e2(1.F, kE2);
    const PGA3D<> e3(1.F, kE3);

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
    px.log();
    printf("a line        : ");
    line.log();
    printf("a plane       : ");
    p.log();
    printf("a rotor       : ");
    rot.log();
    printf("rotated line  : ");
    rotated_line.log();
    printf("rotated point : ");
    rotated_point.log();
    printf("rotated plane : ");
    rotated_plane.log();
    printf("point on plane: ");
    point_on_plane.normalized().log();
    printf("point on torus: ");
    point_on_torus(0.0f, 0.0f).log();
    (e0 - 1.0f).log();
    (1.0f - e0).log();

}

TEST(PGA3DTest, MotorEstimatorTest)
{
    PGA3D<> A = point(0.F, 0.F, 0.F);
    PGA3D<> B = point(1.F, 0.F, 0.F);
    PGA3D<> C = point(0.F,1.F, 0.F);

    PGA3D<> A1 = point(1.F, 1.F, 1.F);
    PGA3D<> B1 = point(1.F, 2.F, 1.F);
    PGA3D<> C1 = point(1.F, 1.F, 2.F);


    // Va = sqrt(Ai/A);
    auto Va = (1.F + (A1*~A)).normalized();
    std::cout << "Va = ";
    Va.log();

    auto Ba = Va * B * ~Va;
    std::cout << "Ba = ";
    Ba.log();

    auto Vb_squared = (A1 & B1)* ~(A1 & Ba);
    std::cout << "Vb_squared = ";
    Vb_squared.log();

    auto Vb = (1 + Vb_squared).normalized();  //Vb = sqrt(Vb_squared);
    std::cout << "Vb = ";
    Vb.log();

    auto Cba = Vb*Va*C*~Va*~Vb;
    std::cout << "Cba = ";
    Cba.log();

    auto Vc_squared = (A1 & B1 & C1) * ~(A1 & B1 & Cba);
    std::cout << "Vc_squared = ";
    Vc_squared.log();

    auto Vc = (1 + Vc_squared).normalized();  //Vc = sqrt(Vc_squared);
    std::cout << "Vc = ";
    Vc.log();

    auto V = Vc * Vb * Va;
    std::cout << "V = ";
    V.log();

    auto Ai = V * A * ~V;
    std::cout << "Ai = ";
    Ai.log();

    auto Bi = V * B * ~V;
    std::cout << "Bi = ";
    Bi.log();

    auto Ci = V * C * ~V;
    std::cout << "Ci = ";
    Ci.log();

}


