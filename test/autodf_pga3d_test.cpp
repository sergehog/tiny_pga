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


#include "../tiny_autodf.h"
#include "../pga3d.h"
#include <gtest/gtest.h>

using namespace tiny_autodf;

using Float = AutoDf<float>;
using APGA = PGA3D<Float>;
using PGA = PGA3D<float>;

TEST(AutoDfPGA3DTest, SimpleTest)
{
    Float::SetType(Float::AutoType::kVariableType);
    Float x = 2.F;
    Float y = 3.F;
    Float z = 4.F;
    Float::SetType(Float::AutoType::kConstType);

    const APGA e0(kE0);
    APGA aa = e0 * x;

    APGA p = point(x, y, z);
    auto p_sq = p*p;
    auto vars1 = p_sq[0].variables();
    EXPECT_EQ(vars1.size(), 3);
    EXPECT_EQ(p_sq[0].value(), -1);

    auto p_sq_dx = p_sq[0].eval();
    EXPECT_EQ(p_sq_dx.derivatives.size(), 3);
}


TEST(AutoDfPGA3DTest, TryTranslatorOptimizationTest)
{
    Float::SetType(Float::AutoType::kConstType);
    APGA p0 = point(Float(0.F), Float(0.F), Float(0.F));
    APGA p1 = point(Float(1.F), Float(1.F), Float(1.F));
    const APGA e01(kE01);
    const APGA e02(kE02);
    const APGA e03(kE03);

    Float::SetType(Float::AutoType::kVariableType);
    Float x = 0.F;
    Float y = 0.F;
    Float z = 0.F;
    Float::SetType(Float::AutoType::kConstType);

    APGA translator = x*e01 + y*e02 + z*e03 + APGA(kScalar);

    auto p_diff = ~(translator*p0*~translator) & ~ p1;
    auto err = p_diff * p_diff;

    auto vars = err[0].variables();
    EXPECT_EQ(vars.size(), 3);

    auto err_dx = err[0].eval();
    EXPECT_EQ(err[0].value(), err_dx.value);
    EXPECT_EQ(err_dx.derivatives.size(), 3);

    size_t i = 0;

    //std::cout << "transl#" << i << ": "; translator.log();
    std::cout << "error#" << i << ": " << err_dx.value << std::endl;


    float err_prev = std::abs(err_dx.value) + 1.F;
    const float learning_rate = -0.01;

    while(std::abs(err_dx.value) < err_prev && i < 20)
    {
        err_prev = std::abs(err_dx.value);
        i ++;
        x.value() -= err_dx.derivatives[x.ID] * learning_rate;
        y.value() -= err_dx.derivatives[y.ID] * learning_rate;
        z.value() -= err_dx.derivatives[z.ID] * learning_rate;

        err_dx = err[0].eval();
        //std::cout << "transl#" << i << ": "; translator.log();
        std::cout << "error#" << i << ": " << err_dx.value << std::endl;
    }


}
