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

//
// Simple automatic differentiation C++/header-only library
// Supposed to be very easy in use, and be able to replace/represent classical scalars types in basic formulas/equations
// Supports partial derivatives (i.e. functions of many variables)
//

#include <cmath>
#include <iostream>
#include <map>
#include <memory>

#ifndef TINY_AUTODF_H_
#define TINY_AUTODF_H_

namespace tiny_autodf
{

template <typename ScalarType = float>
class AutoDf
{
  public:
    struct Evaluation
    {
        ScalarType value;
        std::map<size_t, ScalarType> derivatives;
    };

    enum class AutoType
    {
        kConstType,     // just a const, cannot be changed
        kVariableType,  // input variable, may be assigned/changed
        kSumType,       // result of sum of 2 other AutoDf
        kSubtractType,  // result of subtract of 2 other AutoDf
        kMultType,      // result of multiplication of 2 other AutoDf
        kDivType,       // result of divition of 2 other AutoDf
        kAbsType,
        kMaxType,
        kMinType,
        kSinType,
        kCosType
    };

  private:
    static volatile AutoType default_type_;

    struct CallGraphNode
    {
        size_t ID;
        AutoType type;
        size_t count;
        std::shared_ptr<CallGraphNode> left{};
        std::shared_ptr<CallGraphNode> right{};
        std::shared_ptr<ScalarType> value{};
        std::map<size_t, std::shared_ptr<ScalarType>> variables{};

        CallGraphNode(const size_t id, const AutoType Type, const ScalarType& value_ = static_cast<ScalarType>(0.))
            : ID(id), type(Type)
        {
            value = std::make_shared<ScalarType>(value_);
            count = 1;
        }

        /// Evaluation of call-graph
        /// Returns value & derivatives
        Evaluation eval()
        {
            Evaluation eval_out{};
            if (type == AutoType::kConstType)
            {
                eval_out.value = *value;
                return eval_out;
            }
            else if (type == AutoType::kVariableType)
            {
                eval_out.value = *value;
                eval_out.derivatives[ID] = static_cast<ScalarType>(1.);
                return eval_out;
            }

            const auto eval1 = left->eval();
            const ScalarType v1 = eval1.value;

            if (type == AutoType::kAbsType)
            {
                *value = eval_out.value = std::abs(v1);
                const bool sign_changed = v1 < ScalarType(0.);
                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    const float d1 = eval1.derivatives.at(id);
                    eval_out.derivatives[id] = sign_changed ? -d1 : d1;
                }
                return eval_out;
            }
            else if (type == AutoType::kSinType)
            {
                *value = eval_out.value = std::sin(v1);

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    const float d1 = eval1.derivatives.at(id);
                    eval_out.derivatives[id] = std::cos(eval1.value) * d1;
                }
                return eval_out;
            }
            else if (type == AutoType::kCosType)
            {
                *value = eval_out.value = std::cos(v1);

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    const float d1 = eval1.derivatives.at(id);
                    eval_out.derivatives[id] = -std::sin(eval1.value) * d1;
                }
                return eval_out;
            }

            const auto eval2 = right->eval();
            const ScalarType v2 = eval2.value;

            if (type == AutoType::kSumType)
            {
                *value = eval_out.value = v1 + v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = static_cast<ScalarType>(0.);
                    const auto d1 = eval1.derivatives.find(id);
                    if (d1 != eval1.derivatives.end())
                    {
                        eval_out.derivatives[id] += d1->second;
                    }
                    const auto d2 = eval2.derivatives.find(id);
                    if (d2 != eval2.derivatives.end())
                    {
                        eval_out.derivatives[id] += d2->second;
                    }
                }
            }
            else if (type == AutoType::kSubtractType)
            {
                *value = eval_out.value = v1 - v2;
                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = static_cast<ScalarType>(0.);
                    const auto d1 = eval1.derivatives.find(id);
                    if (d1 != eval1.derivatives.end())
                    {
                        eval_out.derivatives[id] += d1->second;
                    }
                    const auto d2 = eval2.derivatives.find(id);
                    if (d2 != eval2.derivatives.end())
                    {
                        eval_out.derivatives[id] -= d2->second;
                    }
                }
            }
            else if (type == AutoType::kMultType)
            {
                *value = eval_out.value = v1 * v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = static_cast<ScalarType>(0.);
                    const auto d1 = eval1.derivatives.find(id);
                    if (d1 != eval1.derivatives.end())
                    {
                        const ScalarType g1 = d1->second;
                        eval_out.derivatives[id] += v2 * g1;
                    }
                    const auto d2 = eval2.derivatives.find(id);
                    if (d2 != eval2.derivatives.end())
                    {
                        const ScalarType g2 = d2->second;
                        eval_out.derivatives[id] += v1 * g2;
                    }
                }
            }
            else if (type == AutoType::kDivType)
            {
                *value = eval_out.value = v1 / v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    const auto g1 = eval1.derivatives.find(id);
                    const auto g2 = eval2.derivatives.find(id);
                    if (g1 != eval1.derivatives.end())
                    {
                        const ScalarType d1 = g1->second;

                        eval_out.derivatives[id] += d1 * v2 / (v2 * v2);
                    }
                    if (g2 != eval2.derivatives.end())
                    {
                        const ScalarType d2 = g2->second;
                        eval_out.derivatives[id] -= d2 * v1 / (v2 * v2);
                    }
                }
            }
            else if (type == AutoType::kMaxType || type == AutoType::kMinType)
            {
                const bool left_selected = type == AutoType::kMaxType ? v1 >= v2 : v1 <= v2;
                *value = eval_out.value = left_selected ? v1 : v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    const auto d1 = eval1.derivatives.find(id);
                    const auto d2 = eval2.derivatives.find(id);
                    if (d1 != eval1.derivatives.end() && left_selected)
                    {
                        eval_out.derivatives[id] = d1->second;
                    }
                    else
                    {
                        eval_out.derivatives[id] = ScalarType(0.);
                    }
                    if (d2 != eval1.derivatives.end() && !left_selected)
                    {
                        if (d2 != eval1.derivatives.end())
                        {
                            eval_out.derivatives[id] = d2->second;
                        }
                        else
                        {
                            eval_out.derivatives[id] = ScalarType(0.);
                        }
                    }
                }
            }

            return eval_out;
        };
    };

    /// Graph Node belonging to current variable / AutoDf instance
    std::shared_ptr<CallGraphNode> node{};

    /// In order not to multiply zero-valued ConstType nodes, we maintain single such node and re-use
    /// Reduces memory-consumption as well unpredictable id_increment grow
    static const std::shared_ptr<CallGraphNode> zero_node;

    /// Creates AutoDf of specified type, and takes their ownership if needed
    AutoDf(const AutoType type,
           const std::shared_ptr<CallGraphNode>& left,
           const std::shared_ptr<CallGraphNode>& right,
           const ScalarType& value = static_cast<ScalarType>(0.))
    {
        node = std::make_shared<CallGraphNode>((++id_increment), type, value);
        node->left = left;
        node->right = right;
        node->count = 0U;
        if (left)
        {
            node->count += left->count;
            node->variables.insert(left->variables.begin(), left->variables.end());
        }
        if (right)
        {
            node->count += right->count;
            node->variables.insert(right->variables.begin(), right->variables.end());
        }
    }

    AutoDf(const std::shared_ptr<CallGraphNode>& node_in) : node(node_in) {}

    static size_t id_increment;

  public:
    /// Read-only ID value, could be used to distinguish partial derivative of this variable
    size_t ID() const { return node->ID; }

    size_t count() const { return node->count; }

    size_t increment() const { return id_increment; }

    /// Explicitly set what type of AutoDf will be created by-default
    static void StartConstants(const bool need_constant = true)
    {
        default_type_ = (need_constant) ? AutoType::kConstType : AutoType::kVariableType;
    }

    static void StartVariables(const bool need_variable = true)
    {
        default_type_ = (need_variable) ? AutoType::kVariableType : AutoType::kConstType;
    }

    /// Creates kVariableType or kConstType (depending on default_type_) with zero value
    AutoDf()
    {
        if (default_type_ == AutoType::kVariableType)
        {
            node = std::make_shared<CallGraphNode>((++id_increment), AutoType::kVariableType);
            node->variables[node->ID] = node->value;
        }
        else
        {
            node = zero_node;
        }
    }

    /// Creates kVariableType or kConstType (depending on default_type_) with specified value
    AutoDf(const ScalarType& scalar)
    {
        if (default_type_ == AutoType::kConstType && scalar == static_cast<ScalarType>(0))
        {
            node = zero_node;
        }
        else
        {
            node = std::make_shared<CallGraphNode>((++id_increment), default_type_, scalar);
            if (node->type == AutoType::kVariableType)
            {
                node->variables[node->ID] = node->value;
            }
        }
    }

    /// Copy-Constructor, re-uses computation graph node
    AutoDf(const AutoDf<ScalarType>& other) { node = other.node; }

    /// Additional way to create AutoDf of specified type
    AutoDf(const ScalarType& value, const bool is_const)
    {
        const auto type = (is_const ? AutoType::kConstType : AutoType::kVariableType);
        node = std::make_shared<CallGraphNode>((++id_increment), type, value);
        if (node->type == AutoType::kVariableType)
        {
            node->variables[node->ID] = node->value;
        }
    }

    AutoDf& operator=(const ScalarType& scalar)
    {
        if (node == zero_node && scalar != ScalarType(0))
        {
            node = std::make_shared<CallGraphNode>((++id_increment), default_type_, scalar);
            if (default_type_ == AutoType::kVariableType)
            {
                node->variables[node->ID] = node->value;
            }
        }
        else if (node->type == AutoType::kConstType || node->type == AutoType::kVariableType)
        {
            *node->value = scalar;
        }
        return *this;
    }

    AutoDf<ScalarType>& operator=(const AutoDf<ScalarType>& other)
    {
        node = other.node;
        return *this;
    }

    /// List of input "mutable" values, contributing to this AutoDf value
    std::map<size_t, std::shared_ptr<ScalarType>>& variables() const { return node->variables; }

    Evaluation eval() const { return node->eval(); }

    /// Retrieves latest calculated valueabs
    ScalarType& value() const
    {
        if (node->type == AutoType::kVariableType)
        {
            return *node->value;
        }
        else
        {
            // prohibit changing original node value
            static ScalarType fake_value = ScalarType(0);
            fake_value = *node->value;
            return fake_value;
        }
    }

    explicit operator ScalarType() const { return *node->value; }

    AutoDf<ScalarType>& operator+=(const ScalarType value)
    {
        // avoid creating graph nodes, when not really needed
        if (value == static_cast<ScalarType>(0))
        {
            return *this;
        }
        else if (node == zero_node)
        {
            node = std::make_shared<CallGraphNode>((++id_increment), AutoType::kConstType, value);
            return *this;
        }
        else if (node->type == AutoType::kConstType)
        {
            const ScalarType old_node_value = *node->value;
            node = std::make_shared<CallGraphNode>((++id_increment), AutoType::kConstType, old_node_value + value);
            return *this;
        }

        AutoDf<ScalarType> other(value, true);
        return InPlaceOperator(AutoType::kSumType, other.node, *node->value + *other.node->value);
    }

    AutoDf<ScalarType>& operator+=(const AutoDf<ScalarType>& other)
    {
        // avoid creating graph nodes, when not really needed
        if (other.node->type == AutoType::kConstType && *other.node->value == static_cast<ScalarType>(0))
        {
            return *this;
        }
        else if (node->type == AutoType::kConstType && other.node->type == AutoType::kConstType)
        {
            const ScalarType node_value = *node->value;
            node = std::make_shared<CallGraphNode>(
                (++id_increment), AutoType::kConstType, node_value + *other.node->value);
            return *this;
        }

        return InPlaceOperator(AutoType::kSumType, other.node, *node->value + *other.node->value);
    }

    AutoDf<ScalarType>& operator-=(const ScalarType value)
    {
        // avoid creating graph nodes, when not really needed
        if (value == static_cast<ScalarType>(0))
        {
            return *this;
        }
        else if (node->type == AutoType::kConstType)
        {
            const ScalarType node_value = *node->value;
            node = std::make_shared<CallGraphNode>((++id_increment), AutoType::kConstType, node_value - value);
            return *this;
        }

        AutoDf<ScalarType> other(value, true);
        return InPlaceOperator(AutoType::kSubtractType, other.node, *node->value - *other.node->value);
    }

    AutoDf<ScalarType>& operator-=(const AutoDf<ScalarType>& other)
    {
        if (other.node->type == AutoType::kConstType && *other.node->value == static_cast<ScalarType>(0))
        {
            return *this;
        }
        else if (node->type == AutoType::kConstType && other.node->type == AutoType::kConstType)
        {
            const ScalarType node_value = *node->value;
            node = std::make_shared<CallGraphNode>(
                (++id_increment), AutoType::kConstType, node_value - *other.node->value);
            return *this;
        }

        return InPlaceOperator(AutoType::kSubtractType, other.node, *node->value - *other.node->value);
    }

    static AutoDf<ScalarType> abs(const AutoDf<ScalarType>& other)
    {
        return AutoDf<ScalarType>(AutoType::kAbsType, other.node, nullptr, std::abs(*other.node->value));
    }

    static AutoDf<ScalarType> min(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right)
    {
        return AutoDf<ScalarType>(
            AutoType::kMinType, left.node, right.node, std::min(*left.node->value, *right.node->value));
    }

    static AutoDf<ScalarType> max(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right)
    {
        return AutoDf<ScalarType>(
            AutoType::kMaxType, left.node, right.node, std::max(*left.node->value, *right.node->value));
    }

    static AutoDf<ScalarType> sin(const AutoDf<ScalarType>& other)
    {
        return AutoDf<ScalarType>(AutoType::kSinType, other.node, nullptr, std::sin(*other.node->value));
    }

    static AutoDf<ScalarType> cos(const AutoDf<ScalarType>& other)
    {
        return AutoDf<ScalarType>(AutoType::kCosType, other.node, nullptr, std::cos(*other.node->value));
    }

#define DEFINE_OPERATOR(op)                                                     \
    template <typename T>                                                       \
    friend AutoDf<T> operator op(const AutoDf<T>& other, const T scalar_value); \
    template <typename T>                                                       \
    friend AutoDf<T> operator op(const T scalar_value, const AutoDf<T>& other); \
    template <typename T>                                                       \
    friend AutoDf<T> operator op(const AutoDf<T>& left, const AutoDf<T>& right);

    DEFINE_OPERATOR(+);
    DEFINE_OPERATOR(-);
    DEFINE_OPERATOR(*);
    DEFINE_OPERATOR(/);
#undef DEFINE_OPERATOR

    template <typename T>
    friend AutoDf<T> operator-(AutoDf<T> const& other);

  private:
    AutoDf<ScalarType>& InPlaceOperator(const AutoType& type,
                                        const std::shared_ptr<CallGraphNode>& right_node,
                                        const ScalarType new_value)
    {
        auto result = std::make_shared<CallGraphNode>(++id_increment, type, new_value);
        result->left = node;
        result->right = right_node;
        result->count = node->count + right_node->count;
        result->variables.insert(node->variables.begin(), node->variables.end());
        result->variables.insert(right_node->variables.begin(), right_node->variables.end());
        node.swap(result);
        return *this;
    }

    static AutoDf<ScalarType> make_sum(const std::shared_ptr<CallGraphNode>& left,
                                       const std::shared_ptr<CallGraphNode>& right)
    {
        if (left->type == AutoType::kConstType && *left->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(right);
        }
        if (right->type == AutoType::kConstType && *right->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(left);
        }

        const ScalarType value = *left->value + *right->value;
        return AutoDf<ScalarType>(AutoType::kSumType, left, right, value);
    }

    static AutoDf<ScalarType> make_sub(const std::shared_ptr<CallGraphNode>& left,
                                       const std::shared_ptr<CallGraphNode>& right)
    {
        if (right->type == AutoType::kConstType && *right->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(left);
        }

        const ScalarType value = *left->value - *right->value;
        return AutoDf<ScalarType>(AutoType::kSubtractType, left, right, value);
    }

    static AutoDf<ScalarType> make_mult(const std::shared_ptr<CallGraphNode>& left,
                                        const std::shared_ptr<CallGraphNode>& right)
    {
        if (left->type == AutoType::kConstType && *left->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(zero_node);
        }
        if (right->type == AutoType::kConstType && *right->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(zero_node);
        }
        if (left->type == AutoType::kConstType && *left->value == static_cast<ScalarType>(1))
        {
            return AutoDf<ScalarType>(right);
        }
        if (right->type == AutoType::kConstType && *right->value == static_cast<ScalarType>(1))
        {
            return AutoDf<ScalarType>(left);
        }

        const ScalarType value = *left->value * *right->value;
        return AutoDf<ScalarType>(AutoType::kMultType, left, right, value);
    }

    static AutoDf<ScalarType> make_div(const std::shared_ptr<CallGraphNode>& left,
                                       const std::shared_ptr<CallGraphNode>& right)
    {
        if (left->type == AutoType::kConstType && *left->value == static_cast<ScalarType>(0))
        {
            return AutoDf<ScalarType>(zero_node);
        }
        if (right->type == AutoType::kConstType && *right->value == static_cast<ScalarType>(1))
        {
            return AutoDf<ScalarType>(left);
        }

        const ScalarType value = *left->value / *right->value;
        return AutoDf<ScalarType>(AutoType::kDivType, left, right, value);
    }
};

template <typename ScalarType>
AutoDf<ScalarType> operator-(AutoDf<ScalarType> const& other)
{
    if (other.node->type == AutoDf<ScalarType>::AutoType::kConstType)
    {
        return AutoDf<ScalarType>(-other.value(), true);
    }
    else
    {
        return AutoDf<ScalarType>::make_sub(AutoDf<ScalarType>::zero_node, other.node);
    }
}

template <typename ScalarType>
typename AutoDf<ScalarType>::AutoType volatile AutoDf<ScalarType>::default_type_ =
    AutoDf<ScalarType>::AutoType::kVariableType;

#define DEFINE_OPERATOR(OP, func)                                                                   \
    template <typename ScalarType>                                                                  \
    AutoDf<ScalarType> operator OP(const AutoDf<ScalarType>& other, const ScalarType scalar_value)  \
    {                                                                                               \
        AutoDf<ScalarType> scalar(scalar_value, true);                                              \
        return AutoDf<ScalarType>::func(other.node, scalar.node);                                   \
    }                                                                                               \
    template <typename ScalarType>                                                                  \
    AutoDf<ScalarType> operator OP(const ScalarType scalar_value, const AutoDf<ScalarType>& other)  \
    {                                                                                               \
        AutoDf<ScalarType> scalar(scalar_value, true);                                              \
        return AutoDf<ScalarType>::func(scalar.node, other.node);                                   \
    }                                                                                               \
    template <typename ScalarType>                                                                  \
    AutoDf<ScalarType> operator OP(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right) \
    {                                                                                               \
        return AutoDf<ScalarType>::func(left.node, right.node);                                     \
    }

DEFINE_OPERATOR(+, make_sum);
DEFINE_OPERATOR(-, make_sub);
DEFINE_OPERATOR(*, make_mult);
DEFINE_OPERATOR(/, make_div);
#undef DEFINE_OPERATOR

#define INSTANTIATE_AUTODF_TEMPLATE(typename)                                            \
    template <>                                                                          \
    size_t AutoDf<typename>::id_increment = 0U;                                          \
    template <>                                                                          \
    const std::shared_ptr<AutoDf<typename>::CallGraphNode> AutoDf<typename>::zero_node = \
        std::make_shared<AutoDf<typename>::CallGraphNode>(0, AutoType::kConstType, 0.F); \
    AutoDf<float> ____instantiate_template_class;

INSTANTIATE_AUTODF_TEMPLATE(float);

template <typename ScalarType>
typename AutoDf<ScalarType>::Evaluation GradientDescent(const AutoDf<ScalarType>& minimize_error,
                                                        const ScalarType terminate_criteria = ScalarType(1e-6),
                                                        const ScalarType initial_rate = ScalarType(0.01),
                                                        const std::size_t max_iters = 100U)
{
    auto eval = minimize_error.eval();

    if (std::isnan(eval.value) || eval.value < terminate_criteria)
    {
        return eval;
    }

    auto prev_derivatives = eval.derivatives;
    std::map<size_t, float> prev_variables{};
    for (auto x : minimize_error.variables())
    {
        prev_variables[x.first] = *x.second;
    }

    ScalarType prev_error = eval.value + 1.F;
    ScalarType learning_rate = initial_rate;
    size_t iter = 0U;

    while (iter < max_iters && prev_error > terminate_criteria)
    {
        if (iter > 0)
        {
            float dot_value = 0.F;
            float norm_value = 0.F;
            for (auto x : minimize_error.variables())
            {
                const float x_curr = *minimize_error.variables().at(x.first);
                const float x_prev = prev_variables[x.first];
                const float dx_curr = eval.derivatives[x.first];
                const float dx_prev = prev_derivatives[x.first];
                dot_value += (x_curr - x_prev) * (dx_curr - dx_prev);
                norm_value += (dx_curr - dx_prev) * (dx_curr - dx_prev);
            }
            learning_rate = std::abs(dot_value) / norm_value;
        }

        if (std::isnan(learning_rate) || learning_rate < terminate_criteria ||
            std::abs(eval.value) < terminate_criteria)
        {
            return eval;
        }

        std::cout << iter << ": (rate=" << learning_rate << ") F[";

        for (auto dx : eval.derivatives)
        {
            auto x = minimize_error.variables().at(dx.first);
            prev_variables[dx.first] = *x;
            prev_derivatives[dx.first] = dx.second;

            std::cout << *x << ",";
            *x -= dx.second * learning_rate;
        }
        std::cout << "] = " << eval.value << std::endl;
        prev_error = eval.value;
        eval = minimize_error.eval();

        iter++;
    }
    return eval;
}

template <typename ScalarType>
AutoDf<ScalarType> abs(const AutoDf<ScalarType>& other)
{
    return AutoDf<ScalarType>::abs(other);
}

template <typename ScalarType>
AutoDf<ScalarType> min(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right)
{
    return AutoDf<ScalarType>::min(left, right);
}

template <typename ScalarType>
AutoDf<ScalarType> max(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right)
{
    return AutoDf<ScalarType>::max(left, right);
}

template <typename ScalarType>
AutoDf<ScalarType> sin(const AutoDf<ScalarType>& other)
{
    return AutoDf<ScalarType>::sin(other);
}

template <typename ScalarType>
AutoDf<ScalarType> cos(const AutoDf<ScalarType>& other)
{
    return AutoDf<ScalarType>::cos(other);
}

// template <typename ScalarType, std::size_t Size>
// typename AutoDf<ScalarType>::Evaluation LevenbergMarquardt(const std::array<AutoDf<ScalarType>, Size>& A_list,
//                                                          const std::array<ScalarType, Size>& b = {})
//{
//    // list all available variables
//    std::map<size_t, std::shared_ptr<ScalarType>> variables;
//    for(size_t i =0; i<Size; i++)
//    {
//        auto vars = A_list[i].variables();
//        variables.insert(vars.begin(), vars.end());
//    }
//
//}

}  // namespace tiny_autodf

#endif  // TINY_AUTODF_H_
