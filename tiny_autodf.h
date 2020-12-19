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
        kDivType        // result of divition of 2 other AutoDf
    };

  private:
    static AutoType default_type_;

    struct CallGraphNode
    {
        size_t ID;
        AutoType type;
        std::shared_ptr<CallGraphNode> left{};
        std::shared_ptr<CallGraphNode> right{};
        std::shared_ptr<ScalarType> value{};
        std::map<size_t, std::shared_ptr<ScalarType>> variables{};

        CallGraphNode(const size_t id, const AutoType Type, const ScalarType& value_ = static_cast<ScalarType>(0.))
            : ID(id), type(Type)
        {
            value = std::make_shared<ScalarType>(value_);
        }

        Evaluation eval(/*const std::map<size_t, std::shared_ptr<ScalarType>>& values*/)
        {
            Evaluation eval_out{};
            if (type == AutoType::kConstType)
            {
                eval_out.value = *value;
            }
            else if (type == AutoType::kVariableType)
            {
                eval_out.value = *value;
                eval_out.derivatives[ID] = static_cast<ScalarType>(1.);
            }
            else if (type == AutoType::kSumType)
            {
                const auto eval1 = left->eval();
                const auto eval2 = right->eval();
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 + v2;

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
                *value = eval_out.value;
            }
            else if (type == AutoType::kSubtractType)
            {
                const auto eval1 = left->eval();
                const auto eval2 = right->eval();
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 - v2;

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
                *value = eval_out.value;
            }
            else if (type == AutoType::kMultType)
            {
                const auto eval1 = left->eval();
                const auto eval2 = right->eval();
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 * v2;

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
                *value = eval_out.value;
            }
            else if (type == AutoType::kDivType)
            {
                const auto eval1 = left->eval();
                const auto eval2 = right->eval();
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 / v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;

                    const auto d1 = eval1.derivatives.find(id);
                    const auto d2 = eval2.derivatives.find(id);
                    if (d1 != eval1.derivatives.end() && d2 != eval2.derivatives.end())
                    {
                        const ScalarType g1 = d1->second;
                        const ScalarType g2 = d2->second;
                        eval_out.derivatives[id] = (v2 * g1 - v1 * g2) / (g2 * g2);
                    }
                }
                *value = eval_out.value;
            }

            return eval_out;
        };
    };

    // Graph Node belonging to current variable / AutoDf instance
    std::shared_ptr<CallGraphNode> node{};

    /// Creates AutoDf of specified type, and takes their ownership if needed
    AutoDf(const AutoType type,
           const std::shared_ptr<CallGraphNode>& left,
           const std::shared_ptr<CallGraphNode>& right,
           const ScalarType& value = static_cast<ScalarType>(0.))
    {
        node = std::make_shared<CallGraphNode>((++id_increment), type, value);
        node->left = left;
        node->right = right;
    }

    static size_t id_increment;

  public:
    /// Read-only ID value, could be used to distinguish partial derivative of this variable
    size_t ID() { return node->ID; };

    /// Get or Set default type for created AutoDfs
    /// Can be used to distinguish Variables from Consts
    static void SetType(const AutoType& type) { default_type_ = type; }

    /// Empty constructor makes kVariableType with zero value
    AutoDf()
    {
        node = std::make_shared<CallGraphNode>((++id_increment), default_type_);
        if (default_type_ == AutoType::kVariableType)
        {
            node->variables[node->ID] = node->value;
        }
    }

    /// Create kVariableType from ScalarType
    AutoDf(const ScalarType& value)
    {
        node = std::make_shared<CallGraphNode>((++id_increment), default_type_, value);
        if (default_type_ == AutoType::kVariableType)
        {
            node->variables[node->ID] = node->value;
        }
    }

    /// Creating AutoDf as an kCopyTypes
    AutoDf(const AutoDf<ScalarType>& other) { node = other.node; }

    /// Additional way to create kScalarType
    AutoDf(const ScalarType& value, const bool is_const)
    {
        const auto type = (is_const ? AutoType::kConstType : AutoType::kVariableType);
        node = std::make_shared<CallGraphNode>((++id_increment), type, value);
    }

    AutoDf& operator=(const ScalarType& scalar)
    {
        if (node->type == AutoType::kConstType || node->type == AutoType::kVariableType)
        {
            *node->value = scalar;
        }
        else
        {
            // ToDo: throw exception or just ignore?
        }
        return *this;
    }

    AutoDf<ScalarType>& operator=(const AutoDf<ScalarType>& other)
    {
        node = other.node;
        return *this;
    }

    /// List of input "mutable" values, contributing to this AutoDf value
    std::map<size_t, std::shared_ptr<ScalarType>>& variables() { return node->variables; }

    const Evaluation eval() { return node->eval(); }
    //    const Evaluation eval(const std::map<size_t, std::shared_ptr<ScalarType>>& variables)
    //    {
    //        return node->eval(variables);
    //    }

    /// Retrieves latest calculated value
    ScalarType& value() const { return *node->value; }
    explicit operator ScalarType() const { return *node->value; }

    AutoDf<ScalarType>& operator+=(const ScalarType value)
    {
        AutoDf<ScalarType> other(value, true);
        return InPlaceOperator(AutoType::kSumType, other.node, *node->value + *other.node->value);
    }

    AutoDf<ScalarType>& operator+=(const AutoDf<ScalarType>& other)
    {
        return InPlaceOperator(AutoType::kSumType, other.node, *node->value + *other.node->value);
    }

    AutoDf<ScalarType>& operator-=(const ScalarType value)
    {
        AutoDf<ScalarType> other(value, true);
        return InPlaceOperator(AutoType::kSumType, other.node, *node->value - *other.node->value);
    }

    AutoDf<ScalarType>& operator-=(const AutoDf<ScalarType>& other)
    {
        return InPlaceOperator(AutoType::kSumType, other.node, *node->value - *other.node->value);
    }

#define DEFINE_OPERATOR(op)                                                       \
    template <typename T>                                                         \
    friend AutoDf<T>&& operator op(const AutoDf<T>& other, const T scalar_value); \
    template <typename T>                                                         \
    friend AutoDf<T>&& operator op(const T scalar_value, const AutoDf<T>& other); \
    template <typename T>                                                         \
    friend AutoDf<T>&& operator op(const AutoDf<T>& left, const AutoDf<T>& right);

    DEFINE_OPERATOR(+);
    DEFINE_OPERATOR(-);
    DEFINE_OPERATOR(*);
    DEFINE_OPERATOR(/);
#undef DEFINE_OPERATOR

    template <typename T>
    friend AutoDf<T>&& operator-(AutoDf<T> const& other);

  private:
    AutoDf<ScalarType>& InPlaceOperator(const AutoType& type,
                                        const std::shared_ptr<CallGraphNode>& right_node,
                                        const ScalarType new_value)
    {
        auto result = std::make_shared<CallGraphNode>(++id_increment, type, new_value);
        result->left = node;
        result->right = right_node;
        FillVariables(node, right_node, result);
        node.swap(result);
        return *this;
    }

    static AutoDf<ScalarType>&& make_sum(const std::shared_ptr<CallGraphNode>& left,
                                         const std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSumType, left, right);
        *out->node->value = *left->value + *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_sub(const std::shared_ptr<CallGraphNode>& left,
                                         const std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSubtractType, left, right);
        *out->node->value = *left->value - *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_mult(const std::shared_ptr<CallGraphNode>& left,
                                          const std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kMultType, left, right);
        *out->node->value = *left->value * *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_div(const std::shared_ptr<CallGraphNode>& left,
                                         const std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kDivType, left, right);
        *out->node->value = *left->value / *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static void FillVariables(const std::shared_ptr<CallGraphNode>& left,
                              const std::shared_ptr<CallGraphNode>& right,
                              std::shared_ptr<CallGraphNode>& out)
    {
        // copy all values together
        out->variables.insert(left->variables.begin(), left->variables.end());
        out->variables.insert(right->variables.begin(), right->variables.end());
    }
};

template <typename ScalarType>
AutoDf<ScalarType>&& operator-(AutoDf<ScalarType> const& other)
{
    AutoDf<ScalarType> zero(ScalarType(0.0), true);
    return std::move(AutoDf<ScalarType>::make_sub(zero.node, other.node));
}

/// By-default, all created AutoDf considered to be Variables
template <typename ScalarType>
typename AutoDf<ScalarType>::AutoType AutoDf<ScalarType>::default_type_ = AutoDf<ScalarType>::AutoType::kVariableType;

#define DEFINE_OPERATOR(OP, func)                                                                     \
    template <typename ScalarType>                                                                    \
    AutoDf<ScalarType>&& operator OP(const AutoDf<ScalarType>& other, const ScalarType scalar_value)  \
    {                                                                                                 \
        AutoDf<ScalarType> scalar(scalar_value, true);                                                \
        return std::move(AutoDf<ScalarType>::func(other.node, scalar.node));                          \
    }                                                                                                 \
    template <typename ScalarType>                                                                    \
    AutoDf<ScalarType>&& operator OP(const ScalarType scalar_value, const AutoDf<ScalarType>& other)  \
    {                                                                                                 \
        AutoDf<ScalarType> scalar(scalar_value, true);                                                \
        return AutoDf<ScalarType>::func(scalar.node, other.node);                                     \
    }                                                                                                 \
    template <typename ScalarType>                                                                    \
    AutoDf<ScalarType>&& operator OP(const AutoDf<ScalarType>& left, const AutoDf<ScalarType>& right) \
    {                                                                                                 \
        return std::move(AutoDf<ScalarType>::func(left.node, right.node));                            \
    }

DEFINE_OPERATOR(+, make_sum);
DEFINE_OPERATOR(-, make_sub);
DEFINE_OPERATOR(*, make_mult);
DEFINE_OPERATOR(/, make_div);
#undef DEFINE_OPERATOR

template <>
size_t AutoDf<float>::id_increment = 0U;

AutoDf<float> ____instantiate_template_class;

}  // namespace tiny_autodf

#endif  // TINY_AUTODF_H_
