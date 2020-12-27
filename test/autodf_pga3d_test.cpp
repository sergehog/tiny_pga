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
#include "../tiny_autodf.h"
#include <gtest/gtest.h>

using namespace tiny_autodf;
using namespace tiny_pga;
using Float = AutoDf<float>;
using APGA = PGA3D<Float>;
using PGA = PGA3D<float>;

TEST(AutoDfPGA3DTest, SimpleTest)
{
    Float::StartConstants(false);
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float::StartConstants();

    const APGA e0(kE0);
    APGA aa = e0 * x;

    APGA p = point(x, y, z);

    auto p_sq = p[kE013] * p[kE013] + p[kE021] * p[kE021] + p[kE032] * p[kE032];
    auto vars1 = p_sq.variables();
    EXPECT_EQ(vars1.size(), 3);
    EXPECT_EQ(p_sq.value(), 2 * 2 + 3 * 3 + 4 * 4);

    auto p_sq_dx = p_sq.eval();
    EXPECT_EQ(p_sq_dx.derivatives.size(), 3);
}

TEST(AutoDfPGA3DTest, TryTranslatorOptimizationTest)
{
    Float::StartConstants();
    APGA A = point(Float(0.F), Float(0.F), Float(0.F));
    APGA B = point(Float(1.F), Float(0.F), Float(0.F));
    APGA C = point(Float(0.F), Float(1.F), Float(0.F));

    APGA A1 = point(Float(1.F), Float(1.F), Float(1.F));
    APGA B1 = point(Float(1.F), Float(2.F), Float(1.F));
    APGA C1 = point(Float(1.F), Float(1.F), Float(2.F));

    const APGA e12(kE12);
    const APGA e31(kE31);
    const APGA e23(kE23);
    const APGA e01(kE01);
    const APGA e02(kE02);
    const APGA e03(kE03);
    const APGA I(kE0123);

    Float::StartConstants(false);
    Float w = 0.1F;
    Float a = 0.1F;
    Float b = 0.1F;
    Float c = 0.1F;
    Float x = 0.1F;
    Float y = 0.1F;
    Float z = 0.1F;
    Float i = 0.1F;

    Float::StartConstants();

    APGA motor = w * APGA(kScalar) + a * e12 + b * e31 + c * e23 + x * e01 + y * e02 + z * e03 + i * I;

    auto a_diff = motor * A * ~motor - A1;
    auto b_diff = motor * B * ~motor - B1;
    auto c_diff = motor * C * ~motor - C1;

    auto err = a_diff[kE013] * a_diff[kE013] + a_diff[kE021] * a_diff[kE021] + a_diff[kE032] * a_diff[kE032];
    err += b_diff[kE013] * b_diff[kE013] + b_diff[kE021] * b_diff[kE021] + b_diff[kE032] * b_diff[kE032];
    err += c_diff[kE013] * c_diff[kE013] + c_diff[kE021] * c_diff[kE021] + c_diff[kE032] * c_diff[kE032];

    auto vars = err.variables();
    EXPECT_EQ(vars.size(), 8);

    auto eval = err.eval();
    EXPECT_EQ(err.value(), eval.value);
    EXPECT_EQ(eval.derivatives.size(), 8);

    size_t j = 0;

    std::cout << "error#" << j << ": " << eval.value << std::endl;
    std::cout << "motor: [(" << w.value() << "," << a.value() << "," << b.value() << "," << c.value() << "),"
              << x.value() << "," << y.value() << "," << z.value() << "," << i.value() << "]" << std::endl;

    float err_prev = std::abs(eval.value) + 1.F;
    const float learning_rate = 0.02;

    while (std::abs(eval.value) < err_prev && j < 200)
    {
        err_prev = std::abs(eval.value);

        w.value() -= eval.derivatives[w.ID()] * learning_rate;
        a.value() -= eval.derivatives[a.ID()] * learning_rate;
        b.value() -= eval.derivatives[b.ID()] * learning_rate;
        c.value() -= eval.derivatives[c.ID()] * learning_rate;
        x.value() -= eval.derivatives[x.ID()] * learning_rate;
        y.value() -= eval.derivatives[y.ID()] * learning_rate;
        z.value() -= eval.derivatives[z.ID()] * learning_rate;
        i.value() -= eval.derivatives[i.ID()] * learning_rate;
        j++;

        eval = err.eval();

        std::cout << "error#" << j << ": " << eval.value << std::endl;
        std::cout << "motor: [(" << w.value() << "," << a.value() << "," << b.value() << "," << c.value() << "),"
                  << x.value() << "," << y.value() << "," << z.value() << "," << i.value() << "]" << std::endl;
    }

    APGA V = w * APGA(kScalar) + a * e12 + b * e31 + c * e23 + x * e01 + y * e02 + z * e03 + i * I;
    auto Ai = V * A * ~V - A1;
    std::cout << "V*A*~V - A1 = [" << Ai[kE013].value() << "," << Ai[kE021].value() << "," << Ai[kE032].value() << "]"
              << std::endl;

    auto Bi = V * B * ~V - B1;
    std::cout << "V*B*~V - B1 = [" << Bi[kE013].value() << "," << Bi[kE021].value() << "," << Bi[kE032].value() << "]"
              << std::endl;

    auto Ci = V * C * ~V - C1;
    std::cout << "V*C*~V - C1 = [" << Ci[kE013].value() << "," << Ci[kE021].value() << "," << Ci[kE032].value() << "]"
              << std::endl;
}
