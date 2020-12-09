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

#ifndef TINY_AUTODIFF_H_
#define TINY_AUTODIFF_H_

namespace tiny_autodiff
{

template <typename ScalarType = float>
class AutoDf
{
  private:
    enum class AutoType
    {
        kScalarType,    // just a const, cannot be changed
        kVariableType,  // input variable, may be assigned/changed
        kCopyType,      // result of copy/lvalue assignment/etc from other AutoDf
        kSumType,       // result of sum of 2 other AutoDf
        kSubtractType,  // result of subtract of 2 other AutoDf
        kMultType,      // result of multiplication of 2 other AutoDf
        kDivType        // result of divition of 2 other AutoDf
    };

    const AutoType type_;
    static size_t id_increment;

    AutoDf<ScalarType>* left_ = nullptr;
    AutoDf<ScalarType>* right_ = nullptr;
    const bool own_left_ = false;
    const bool own_right_ = false;

    std::shared_ptr<ScalarType> value_{};
    std::map<size_t, std::shared_ptr<ScalarType>> variables_{};

    /// Creates AutoDf of specified type, and takes their ownership if needed
    AutoDf(AutoType type,
           AutoDf<ScalarType>* left = nullptr,
           AutoDf<ScalarType>* right = nullptr,
           const bool own_left = false,
           const bool own_right = false)
        : ID(++id_increment), type_(type), own_left_(left != nullptr && own_left), own_right_(right != nullptr && own_right)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        left_ = left;
        right_ = right;
    }

  public:

    /// Read-only ID value
    /// Could be used to distinguish it's derivative
    const size_t ID;

    /// Holds estimated value and list of derivatives
    struct Evaluation
    {
        ScalarType value;
        std::map<size_t, ScalarType> derivatives;
    };


    /// Empty constructor makes kVariableType with zero value
    AutoDf() : ID(++id_increment), type_(AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = 0.;
        variables_[ID] = value_;
    }

    /// Creating AutoDf as an kVariableType or kScalarType
    explicit AutoDf(const ScalarType&& value)
        : ID(++id_increment), type_(AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = value;
        variables_[ID] = value_;
    }

    /// Creating AutoDf as an kVariableType or kScalarType
    explicit AutoDf(const ScalarType& value)
        : ID(++id_increment), type_(AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = value;
        variables_[ID] = value_;
    }

//    AutoDf<ScalarType>&& operator=(const ScalarType value)
//    {
//        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(value);
//        return std::move(*out);
//    }
//
//    AutoDf<ScalarType>&& operator=(const ScalarType& value)
//    {
//        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(value);
//        return std::move(*out);
//    }

    /// Creating AutoDf as an kVariableType or kScalarType
    AutoDf(const ScalarType& value, bool is_const)
        : ID(++id_increment), type_(is_const ? AutoType::kScalarType : AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = value;

        if (type_ == AutoType::kVariableType)
        {
            variables_[ID] = value_;
        }
    }

//    /// Creating AutoDf as an kCopyValue
//    explicit AutoDf(AutoDf<ScalarType>& other) : ID{++id_increment}, type_(AutoType::kCopyType), left_(&other)
//    {
//        value_ = std::shared_ptr<ScalarType>(new ScalarType);
//        *value_ = *left_->value_;
//        std::copy(other.variables_.begin(), other.variables_.end(), variables_.begin());
//    }

    explicit AutoDf(AutoDf<ScalarType>&& other)
        : ID{other.ID},
          type_(other.type_),
          left_(other.left_),
          right_(other.right_),
          own_left_(other.own_left_),
          own_right_(other.own_right_)
    {
        value_ = std::move(other.value_);
        variables_ = std::move(other.variables_);
    }

    /// List of input "mutable" values, contributing to this AutoDf value
    std::map<size_t, std::shared_ptr<ScalarType>>& variables() { return variables_; }

    Evaluation eval(std::map<size_t, std::shared_ptr<ScalarType>> values)
    {
        Evaluation eval_out{};
        if (type_ == AutoType::kScalarType)
        {
            eval_out.value = *value_;
        }
        else if (type_ == AutoType::kVariableType)
        {
            eval_out.value = *value_;
            eval_out.derivatives[ID] = 1.;
        }
        else if (type_ == AutoType::kCopyType)
        {
            eval_out = left_->eval(values);
            *value_ = eval_out.value;
        }
        else if (type_ == AutoType::kSumType)
        {
            const auto eval1 = left_->eval(values);
            const auto eval2 = right_->eval(values);
            const ScalarType v1 = eval1.value;
            const ScalarType v2 = eval2.value;
            eval_out.value = v1 + v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
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
            *value_ = eval_out.value;
        }
        else if (type_ == AutoType::kSubtractType)
        {
            const auto eval1 = left_->eval(values);
            const auto eval2 = right_->eval(values);
            const ScalarType v1 = eval1.value;
            const ScalarType v2 = eval2.value;
            eval_out.value = v1 - v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
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
            *value_ = eval_out.value;
        }
        else if (type_ == AutoType::kMultType)
        {
            const auto eval1 = left_->eval(values);
            const auto eval2 = right_->eval(values);
            const ScalarType v1 = eval1.value;
            const ScalarType v2 = eval2.value;
            eval_out.value = v1 * v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
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
            *value_ = eval_out.value;
        }
        else if (type_ == AutoType::kDivType)
        {
            const auto eval1 = left_->eval(values);
            const auto eval2 = right_->eval(values);
            const ScalarType v1 = eval1.value;
            const ScalarType v2 = eval2.value;
            eval_out.value = v1 / v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
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
            *value_ = eval_out.value;
        }

        return eval_out;
    };

    Evaluation eval() { return eval(variables_); }

    ScalarType& value() const { return *value_; }

    AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>& other)
    {
        return make_sum(this, &other);
    }

    AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>&& other)
    {
        AutoDf<float>* copy = new AutoDf<float>(std::move(other));
        return make_sum(this, *copy, false, true);
    }

    AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>& other)
    {
        return make_sub(this, &other);
    }

    AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>&& other)
    {
        AutoDf<float>* copy = new AutoDf<float>(std::move(other));
        return make_sub(this, copy, false, true);
    }

    AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>& other)
    {
        return make_mult(this, &other);
    }

    AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>&& other)
    {
        AutoDf<float>* copy = new AutoDf<float>(std::move(other));
        return make_mult(this, copy, false, true);
    }

    AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>& other)
    {
        return make_div(this, &other);
    }

    AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>&& other)
    {
        AutoDf<float>* copy = new AutoDf<float>(std::move(other));
        return make_div(this, copy, false, true);
    }

    ~AutoDf()
    {
        if (own_left_ && left_ != nullptr)
        {
            delete left_;
            left_ = nullptr;
        }
        if (own_right_ && right_ != nullptr)
        {
            delete right_;
            right_ = nullptr;
        }
        variables_.clear();
        value_.reset();
    }

    /// addition from the right
    friend AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType>&& operator+(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// addition from the left
    friend AutoDf<ScalarType>&& operator+(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    //friend AutoDf<ScalarType>&& operator+(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    /// subtract from the right
    friend AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    //friend AutoDf<ScalarType>&& operator-(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// subtract from the left
    friend AutoDf<ScalarType>&& operator-(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    //friend AutoDf<ScalarType>&& operator-(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    /// multiplication from the left
    friend AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    //friend AutoDf<ScalarType>&& operator*(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    /// multiplication from the right
    friend AutoDf<ScalarType>&& operator*(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    //friend AutoDf<ScalarType>&& operator*(const ScalarType scalar_value, AutoDf<ScalarType>&& other);

    friend AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    //friend AutoDf<ScalarType>&& operator/(AutoDf<ScalarType>&& other, const ScalarType scalar_value);

    friend AutoDf<ScalarType>&& operator/(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    //friend AutoDf<ScalarType>&& operator/(const ScalarType scalar_value, AutoDf<ScalarType>&& other);


  private:
    static AutoDf<ScalarType>&& make_sum(AutoDf<ScalarType>* left,
                                       AutoDf<ScalarType>* right,
                                       const bool own_left = false,
                                       const bool own_right = false)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSumType, left, right, own_left, own_right);
        *out->value_ = *left->value_ + *right->value_;
        FillVariables(left, right, out);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_sub(AutoDf<ScalarType>* left,
                                       AutoDf<ScalarType>* right,
                                       const bool own_left = false,
                                       const bool own_right = false)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kSubtractType, left, right, own_left, own_right);
        *out->value_ = *left->value_ - *right->value_;
        FillVariables(left, right, out);
        return std::move(*out);
    }

    static AutoDf<ScalarType>&& make_mult(AutoDf<ScalarType>* left,
                                        AutoDf<ScalarType>* right,
                                        const bool own_left = false,
                                        const bool own_right = false)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kMultType, left, right, own_left, own_right);
        *out->value_ = *left->value_ * *right->value_;
        FillVariables(left, right, out);
        return std::move(*out);
    }


    static AutoDf<ScalarType>&& make_div(AutoDf<ScalarType>* left,
                                       AutoDf<ScalarType>* right,
                                       const bool own_left = false,
                                       const bool own_right = false)
    {
        AutoDf<ScalarType>* out = new AutoDf<ScalarType>(AutoType::kDivType, left, right, own_left, own_right);
        *out->value_ = *left->value_ / *right->value_;
        FillVariables(left, right, out);
        return std::move(*out);
    }

    static void FillVariables(AutoDf<ScalarType>* left,
                              AutoDf<ScalarType>* right,
                              AutoDf<ScalarType>* out)
    {
        // copy all values together
        for (auto v = left->variables_.begin(); v != left->variables_.end(); v++)
        {
            out->variables_[v->first] = v->second;
        }
        for (auto v = right->variables_.begin(); v != right->variables_.end(); v++)
        {
            out->variables_[v->first] = v->second;
        }
    }
};

AutoDf<float>&& operator+(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(&other, scalar, false, true));
}

AutoDf<float>&& operator+(AutoDf<float>&& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_sum(&other, scalar, false, true));
}

AutoDf<float>&& operator+(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_sum( scalar, &other, true));
}

AutoDf<float>&& operator-(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_sub(&other, scalar, false, true));
}

AutoDf<float>&& operator-(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_sub(scalar, &other, true));
}

AutoDf<float>&& operator*(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_mult(&other, scalar, false, true));
}

AutoDf<float>&& operator*(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_mult(scalar, &other, true, false));
}

AutoDf<float>&& operator/(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_div(&other, scalar, false, true));
}

AutoDf<float>&& operator/(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    return std::move(AutoDf<float>::make_div(scalar, &other, true, false));
}


template <>
size_t AutoDf<float>::id_increment = 0U;


AutoDf<float> ____instantiate_template_class;

}  // namespace tiny_autodiff

#endif  // TINY_AUTODIFF_H_
