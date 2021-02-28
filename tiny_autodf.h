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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

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
        std::unordered_map<size_t, ScalarType> derivatives;
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
        std::unordered_map<size_t, std::shared_ptr<ScalarType>> variables{};

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
    std::unordered_map<size_t, std::shared_ptr<ScalarType>>& variables() const { return node->variables; }

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

/// MACRO to helps define operators of 2 arguments
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

/// Definition of math operators
DEFINE_OPERATOR(+, make_sum);
DEFINE_OPERATOR(-, make_sub);
DEFINE_OPERATOR(*, make_mult);
DEFINE_OPERATOR(/, make_div);
#undef DEFINE_OPERATOR

/// MACRO to helps define functions of 2 arguments
#define DEFINE_FUNCTION(func)                                                                \
    template <typename ScalarType>                                                           \
    AutoDf<ScalarType> func(const AutoDf<ScalarType>& other, const ScalarType scalar_value)  \
    {                                                                                        \
        AutoDf<ScalarType> scalar(scalar_value, true);                                       \
        return AutoDf<ScalarType>::func(other, scalar);                                      \
    }                                                                                        \
    template <typename ScalarType>                                                           \
    AutoDf<ScalarType> func(const ScalarType scalar_value, const AutoDf<ScalarType>& other)  \
    {                                                                                        \
        AutoDf<ScalarType> scalar(scalar_value, true);                                       \
        return AutoDf<ScalarType>::func(scalar, other);                                      \
    }                                                                                        \
    template <typename ScalarType>                                                           \
    AutoDf<ScalarType> func(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right) \
    {                                                                                        \
        return AutoDf<ScalarType>::func(left, right);                                        \
    }

/// Definition of some functions
DEFINE_FUNCTION(min);
DEFINE_FUNCTION(max);
#undef DEFINE_FUNCTION

/// Functions of 1 argument do not require macros
template <typename ScalarType>
AutoDf<ScalarType> abs(const AutoDf<ScalarType>& other)
{
    return AutoDf<ScalarType>::abs(other);
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

#define INSTANTIATE_AUTODF_TEMPLATE(typename)                                            \
    template <>                                                                          \
    size_t AutoDf<typename>::id_increment = 0U;                                          \
    template <>                                                                          \
    const std::shared_ptr<AutoDf<typename>::CallGraphNode> AutoDf<typename>::zero_node = \
        std::make_shared<AutoDf<typename>::CallGraphNode>(0, AutoType::kConstType, 0.F); \
    AutoDf<typename> ____instantiate_AutoDf_##typename;

INSTANTIATE_AUTODF_TEMPLATE(float);

/// In order to instantiate AutoDf<double> class one need either uncomment next line, or write one in their code
// INSTANTIATE_AUTODF_TEMPLATE(double);

template <typename ScalarType>
struct TerminationCriteria
{
    ScalarType expression_less_than = NAN;
    ScalarType step_less_than = 1e-6;
    ScalarType diff_less_than = 1e-8;
};

/// Gradient Descent minimization algorithm
template <typename ScalarType>
typename AutoDf<ScalarType>::Evaluation GradientDescent(
    const AutoDf<ScalarType>& minimize_expression,
    const TerminationCriteria<ScalarType> termination_criteria = {NAN, 1e-6, 1e-8},
    const ScalarType initial_step = ScalarType(0.01),
    const std::size_t max_iterations = 100U)
{
    auto prev_error = minimize_expression.eval();
    const bool termination_criteria_met = (!std::isnan(termination_criteria.expression_less_than)) &&
                                          (prev_error.value < termination_criteria.expression_less_than);

    if (std::isnan(prev_error.value) || termination_criteria_met)
    {
        return prev_error;
    }

    std::cout << "0: (initial) F[.]=" << prev_error.value << std::endl;

    // store previous variable values
    std::unordered_map<size_t, float> prev_values{};
    for (auto value_pair : minimize_expression.variables())
    {
        prev_values[value_pair.first] = *value_pair.second;
    }

    ScalarType current_step = initial_step;

    for (size_t iter = 1U; iter < max_iterations; iter++)
    {
        std::cout << iter << ": (rate=" << current_step << ") F[";

        for (auto value_pair : minimize_expression.variables())
        {
            const ScalarType derivative = prev_error.derivatives[value_pair.first];
            *value_pair.second = prev_values[value_pair.first] - derivative * current_step;
            std::cout << *value_pair.second << ",";
        }

        auto current_error = minimize_expression.eval();
        std::cout << "] = " << current_error.value << std::endl;

        const bool expression_criteria_met = (!std::isnan(termination_criteria.expression_less_than)) &&
                                             (current_error.value < termination_criteria.expression_less_than);

        ScalarType error_diff = prev_error.value - current_error.value;
        const bool diff_criteria_met =
            (!std::isnan(termination_criteria.diff_less_than)) && (error_diff < termination_criteria.diff_less_than);

        if (expression_criteria_met || diff_criteria_met)
        {
            return current_error;
        }

        // if error get worse -> "re-play" previous iteration with smaller step
        if (current_error.value > prev_error.value)
        {
            // reduce step
            current_step /= 2.;
        }
        else
        {

            float dot_value = 0.F, norm_value = 0.F;

            for (auto value_pair : minimize_expression.variables())
            {
                const ScalarType x_prev = prev_values[value_pair.first];
                const ScalarType dx_prev = prev_error.derivatives[value_pair.first];

                ScalarType x_curr = *value_pair.second;
                const ScalarType dx_curr = current_error.derivatives[value_pair.first];

                dot_value += (x_curr - x_prev) * (dx_curr - dx_prev);
                norm_value += (dx_curr - dx_prev) * (dx_curr - dx_prev);

                // update previous values
                prev_values[value_pair.first] = x_curr;
            }

            current_step = std::abs(dot_value) / norm_value;
            prev_error = current_error;
        }

        const bool step_criteria_met =
            (!std::isnan(termination_criteria.step_less_than)) && (current_step < termination_criteria.step_less_than);

        if (step_criteria_met)
        {
            // if current_step is getting too small -> restore previous values and return previuous error
            for (auto value_pair : minimize_expression.variables())
            {
                *value_pair.second = prev_values[value_pair.first];
            }
            return prev_error;
        }
    }

    // return best known error
    for (auto value_pair : minimize_expression.variables())
    {
        *value_pair.second = prev_values[value_pair.first];
    }
    return prev_error;
}
//
// template <typename ScalarType, std::size_t Size>
// typename AutoDf<ScalarType>::Evaluation LevenbergMarquardt(const std::array<AutoDf<ScalarType>, Size>& A_list,
//                                                          const std::array<ScalarType, Size>& b = {})
//{
//    // list all available variables
//    std::map<size_t, std::shared_ptr<ScalarType>> variables {};
//    std::for_each(A_list.begin(), A_list.end(), [&](auto& a){variables.insert(a.variables().begin());});
//
//}

}  // namespace tiny_autodf

#endif  // TINY_AUTODF_H_
