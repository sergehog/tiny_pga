//
// Simple automatic differentiation C++/header-only library
// Supposed to be very easy in use, and be able to replace/represent classical scalars types in basic formulas/equations
// Supports partial derivatives (i.e. functions of many variables)
//
// Copyright Sergey Smirnov  / Seregium Oy 2020
// Email: sergei.smirnov@gmail.com
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

  private:
    enum class AutoType
    {
        kScalarType,    // just a const, cannot be changed
        kVariableType,  // input variable, may be assigned/changed
        // kCopyType,      // result of copy/lvalue assignment/etc from other AutoDf
        kSumType,       // result of sum of 2 other AutoDf
        kSubtractType,  // result of subtract of 2 other AutoDf
        kMultType,      // result of multiplication of 2 other AutoDf
        kDivType        // result of divition of 2 other AutoDf
    };

    struct CallGraphNode
    {
        size_t ID;
        AutoType type;
        std::shared_ptr<CallGraphNode> left{};
        std::shared_ptr<CallGraphNode> right{};
        std::shared_ptr<ScalarType> value{};
        std::map<size_t, std::shared_ptr<ScalarType>> variables{};

        CallGraphNode(const size_t id, AutoType Type, const ScalarType& value_ = 0) : ID(id), type(Type)
        {
            value = std::make_shared<ScalarType>(value_);
        }

        Evaluation eval(const std::map<size_t, std::shared_ptr<ScalarType>>& values)
        {
            Evaluation eval_out{};
            if (type == AutoType::kScalarType)
            {
                eval_out.value = *value;
            }
            else if (type == AutoType::kVariableType)
            {
                eval_out.value = *value;
                eval_out.derivatives[ID] = 1.;
            }
            //        else if (node->type == AutoType::kCopyType)
            //        {
            //            eval_out = node->left->eval(values);
            //            *node->value = eval_out.value;
            //        }
            else if (type == AutoType::kSumType)
            {
                const auto eval1 = left->eval(values);
                const auto eval2 = right->eval(values);
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 + v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = 0.F;
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
                const auto eval1 = left->eval(values);
                const auto eval2 = right->eval(values);
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 - v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = 0.F;
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
                const auto eval1 = left->eval(values);
                const auto eval2 = right->eval(values);
                const ScalarType v1 = eval1.value;
                const ScalarType v2 = eval2.value;
                eval_out.value = v1 * v2;

                for (auto vi = variables.begin(); vi != variables.end(); vi++)
                {
                    const size_t id = vi->first;
                    eval_out.derivatives[id] = 0.F;
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
                const auto eval1 = left->eval(values);
                const auto eval2 = right->eval(values);
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
           std::shared_ptr<CallGraphNode>& left,
           std::shared_ptr<CallGraphNode>& right,
           const ScalarType& value = 0)
        : ID(++id_increment)
    {
        node = std::make_shared<CallGraphNode>(ID, type, value);
        node->left = left;
        node->right = right;
    }

    static size_t id_increment;

  public:

    /// Read-only ID value, could be used to distinguish its derivative
    const size_t ID;

    /// Empty constructor makes kVariableType with zero value
    AutoDf() : ID(++id_increment)
    {
        node = std::make_shared<CallGraphNode>(ID, AutoType::kVariableType, 0.);
        node->variables[node->ID] = node->value;
    }

    /// Create kVariableType from ScalarType
    AutoDf(const ScalarType& value) : ID(++id_increment)
    {
        node = std::make_shared<CallGraphNode>(ID, AutoType::kVariableType, value);
        node->variables[node->ID] = node->value;
    }

//    /// Moving Ctor from rvalue -> takes ownership of everything
//    AutoDf(AutoDf<ScalarType>&& other) : ID(other.ID) { node = std::move(other.node); }

    /// Creating AutoDf as an kCopyTypes
     AutoDf(const AutoDf<ScalarType>& other) : ID(++id_increment) { node = other.node; }

    /// Additional way to create kScalarType
    AutoDf(const ScalarType& value, const bool is_const) : ID(++id_increment)
    {
        const auto type = (is_const ? AutoType::kScalarType : AutoType::kVariableType);
        node = std::make_shared<CallGraphNode>(ID, type, value);
    }

    /// List of input "mutable" values, contributing to this AutoDf value
    std::map<size_t, std::shared_ptr<ScalarType>>& variables() { return node->variables; }

    const Evaluation& eval() { return node->eval(node->variables); }
    const Evaluation eval(const std::map<size_t, std::shared_ptr<ScalarType>>& variables)
    {
        return node->eval(variables);
    }

    /// Retrieves latest calculated value
    ScalarType& value() const { return *node->value; }
    operator ScalarType() const { return *node->value; }

    AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>& other) { return make_sum(node, other.node); }
    AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>&& other) { return make_sum(node, other.node); }
    AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>& other) { return make_sub(node, other.node); }
    AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>&& other) { return make_sub(node, other.node); }
    AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>& other) { return make_mult(node, other.node); }
    AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>&& other) { return make_mult(node, other.node); }
    AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>& other) { return make_div(node, other.node); }
    AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>&& other) { return make_div(node, other.node); }

    //    AutoDf<ScalarType>& operator+=(AutoDf<ScalarType>& other)
    //    {
    //        return make_sum(this, &other);
    //    }

    /// addition from the right
    friend AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// addition from the left
    friend AutoDf<ScalarType>&& operator+(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType>&& operator+(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    /// subtract from the right
    friend AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// subtract from the left
    friend AutoDf<ScalarType>&& operator-(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType>&& operator-(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    /// multiplication from the left
    friend AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// multiplication from the right
    friend AutoDf<ScalarType>&& operator*(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType>&& operator*(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    friend AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    friend AutoDf<ScalarType>&& operator/(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType>&& operator/(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

  private:
    static AutoDf<ScalarType>&& make_sum(std::shared_ptr<CallGraphNode>& left, std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSumType, left, right);
        *out->node->value = *left->value + *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_sub(std::shared_ptr<CallGraphNode>& left, std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSubtractType, left, right);
        *out->node->value = *left->value - *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_mult(std::shared_ptr<CallGraphNode>& left, std::shared_ptr<CallGraphNode>& right)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kMultType, left, right);
        *out->node->value = *left->value * *right->value;
        FillVariables(left, right, out->node);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_div(std::shared_ptr<CallGraphNode>& left, std::shared_ptr<CallGraphNode>& right)
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
        for (auto v = left->variables.begin(); v != left->variables.end(); v++)
        {
            out->variables[v->first] = v->second;
        }
        for (auto v = right->variables.begin(); v != right->variables.end(); v++)
        {
            out->variables[v->first] = v->second;
        }
    }
};

AutoDf<float>&& operator+(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(other.node, scalar.node));
}

AutoDf<float>&& operator+(AutoDf<float>&& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(other.node, scalar.node));
}

AutoDf<float>&& operator+(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(scalar.node, other.node));
}

AutoDf<float>&& operator+(const float scalar_value, AutoDf<float>&& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(scalar.node, other.node));
}

AutoDf<float>&& operator-(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sub(other.node, scalar.node));
}
AutoDf<float>&& operator-(AutoDf<float>&& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sub(other.node, scalar.node));
}

AutoDf<float>&& operator-(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_sub(scalar.node, other.node));
}

AutoDf<float>&& operator-(const float scalar_value, AutoDf<float>&& other)
{
    AutoDf<float> scalar(scalar_value, true);
    AutoDf<float>* other_copy = new AutoDf<float>(std::move(other));
    return std::move(AutoDf<float>::make_sub(scalar.node, other.node));
}

AutoDf<float>&& operator*(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_mult(other.node, scalar.node));
}

AutoDf<float>&& operator*(AutoDf<float>&& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    AutoDf<float>* other_copy = new AutoDf<float>(std::move(other));
    return std::move(AutoDf<float>::make_mult(other.node, scalar.node));
}

AutoDf<float>&& operator*(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_mult(scalar.node, other.node));
}

AutoDf<float>&& operator*(const float scalar_value, AutoDf<float>&& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_mult(scalar.node, other.node));
}

AutoDf<float>&& operator/(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_div(other.node, scalar.node));
}

AutoDf<float>&& operator/(AutoDf<float>&& other, const float scalar_value)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_div(other.node, scalar.node));
}

AutoDf<float>&& operator/(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_div(scalar.node, other.node));
}

AutoDf<float>&& operator/(const float scalar_value, AutoDf<float>&& other)
{
    AutoDf<float> scalar(scalar_value, true);
    return std::move(AutoDf<float>::make_div(scalar.node, other.node));
}

template <>
size_t AutoDf<float>::id_increment = 0U;

AutoDf<float> ____instantiate_template_class;

}  // namespace tiny_autodf

#endif  // TINY_AUTODF_H_
