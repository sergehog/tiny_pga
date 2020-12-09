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
        kCopyType,      // result of copy/assignment from other AutoDf
        kSumType,       // result of sum of 2 other AutoDf
        kSubtractType,  // result of subtract of 2 other AutoDf
        kMultType,      // result of multiplication of 2 other AutoDf
        kDivType        // result of divition of 2 other AutoDf
    };

    const AutoType type_;
    const size_t id_;
    static size_t id_increment;

    std::shared_ptr<ScalarType> value_{};

    AutoDf<ScalarType>* first_ = nullptr;
    AutoDf<ScalarType>* second_ = nullptr;
    const bool own_first = false;
    const bool own_second = false;

    std::map<size_t, std::shared_ptr<ScalarType>> variables_{};

    /// Creating AutoDf as an kInputValue or kConstValue
    /// take ownership of provided pointers
    AutoDf(AutoType type, AutoDf<ScalarType>* first = nullptr, AutoDf<ScalarType>* second = nullptr)
        : id_(0U), type_(type), own_first(first != nullptr), own_second(second != nullptr)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        first_ = first;
        second_ = second;
    }

  public:
    /// Holds estimated value and list of gradiens
    using EvaluationType = std::pair<ScalarType, std::map<size_t, ScalarType>>;

    /// Empty constructor makes kInputValue with zero value
    AutoDf() : id_(++id_increment), type_(AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = 0.;
        variables_[id_] = value_;
    }

    /// Creating AutoDf as an kVariableType or kScalarType
    AutoDf(const ScalarType& const_value, bool is_const = false)
        : id_(is_const ? 0U : ++id_increment), type_(is_const ? AutoType::kScalarType : AutoType::kVariableType)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        *value_ = const_value;

        if (type_ == AutoType::kVariableType)
        {
            variables_[id_] = value_;
        }
    }

    /// Creating AutoDf as an kCopyValue
    AutoDf(AutoDf<ScalarType>& other) : id_{}, type_(AutoType::kCopyType), first_(&other)
    {
        value_ = std::shared_ptr<ScalarType>(new ScalarType);
        std::copy(other.variables_.begin(), other.variables_.end(), variables_.begin());
    }

    AutoDf(AutoDf<ScalarType>&& other)
        : id_{other.id_},
          type_(other.type_),
          first_(other.first_),
          second_(other.second_),
          own_first(other.own_first),
          own_second(other.own_second)
    {
        value_ = std::move(other.value_);
        variables_ = std::move(other.variables_);
    }

    /// List of input "mutable" values, contributing to this AutoDf value
    std::map<size_t, std::shared_ptr<ScalarType>>& variables() { return variables_; }

    EvaluationType eval(std::map<size_t, std::shared_ptr<ScalarType>> values)
    {
        EvaluationType eval_out{};
        if (type_ == AutoType::kScalarType)
        {
            eval_out.first = *value_;
        }
        else if (type_ == AutoType::kVariableType)
        {
            eval_out.first = *value_;
            eval_out.second[id_] = 1.;
        }
        else if (type_ == AutoType::kCopyType)
        {
            eval_out = first_->eval(values);
            *value_ = eval_out.first;
        }
        else if (type_ == AutoType::kSumType)
        {
            const auto eval1 = first_->eval(values);
            const auto eval2 = second_->eval(values);
            const ScalarType v1 = eval1.first;
            const ScalarType v2 = eval2.first;
            eval_out.first = v1 + v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
            {
                const size_t id = vi->first;
                eval_out.second[id] = 0.F;
                const auto d1 = eval1.second.find(id);
                if (d1 != eval1.second.end())
                {
                    eval_out.second[id] += d1->second;
                }
                const auto d2 = eval2.second.find(id);
                if (d2 != eval2.second.end())
                {
                    eval_out.second[id] += d2->second;
                }
            }
            *value_ = eval_out.first;
        }
        else if (type_ == AutoType::kSubtractType)
        {
            const auto eval1 = first_->eval(values);
            const auto eval2 = second_->eval(values);
            const ScalarType v1 = eval1.first;
            const ScalarType v2 = eval2.first;
            eval_out.first = v1 - v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
            {
                const size_t id = vi->first;
                eval_out.second[id] = 0.F;
                const auto d1 = eval1.second.find(id);
                if (d1 != eval1.second.end())
                {
                    eval_out.second[id] += d1->second;
                }
                const auto d2 = eval2.second.find(id);
                if (d2 != eval2.second.end())
                {
                    eval_out.second[id] -= d2->second;
                }
            }
            *value_ = eval_out.first;
        }
        else if (type_ == AutoType::kMultType)
        {
            const auto eval1 = first_->eval(values);
            const auto eval2 = second_->eval(values);
            const ScalarType v1 = eval1.first;
            const ScalarType v2 = eval2.first;
            eval_out.first = v1 * v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
            {
                const size_t id = vi->first;
                eval_out.second[id] = 0.F;
                const auto d1 = eval1.second.find(id);
                if (d1 != eval1.second.end())
                {
                    const ScalarType g1 = d1->second;
                    eval_out.second[id] += v2 * g1;
                }
                const auto d2 = eval2.second.find(id);
                if (d2 != eval2.second.end())
                {
                    const ScalarType g2 = d2->second;
                    eval_out.second[id] += v1 * g2;
                }
            }
            *value_ = eval_out.first;
        }
        else if (type_ == AutoType::kDivType)
        {
            const auto eval1 = first_->eval(values);
            const auto eval2 = second_->eval(values);
            const ScalarType v1 = eval1.first;
            const ScalarType v2 = eval2.first;
            eval_out.first = v1 / v2;

            for (auto vi = variables_.begin(); vi != variables_.end(); vi++)
            {
                const size_t id = vi->first;

                const auto d1 = eval1.second.find(id);
                const auto d2 = eval2.second.find(id);
                if (d1 != eval1.second.end() && d2 != eval2.second.end())
                {
                    const ScalarType g1 = d1->second;
                    const ScalarType g2 = d2->second;
                    eval_out.second[id] = (v2 * g1 - v1 * g2) / (g2 * g2);
                }
            }
            *value_ = eval_out.first;
        }

        return eval_out;
    };

    EvaluationType eval() { return eval(variables_); }

    ScalarType value() const { return *value_; }

    AutoDf<ScalarType> operator+(AutoDf<ScalarType>& other)
    {
        AutoDf<ScalarType> mult(AutoType::kSumType);
        mult.first_ = this;
        mult.second_ = &other;
        *mult.value_ = *value_ + *other.value_;

        // copy all values together
        for (auto v = variables_.begin(); v != variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        // overwrites if some of them the same
        for (auto v = other.variables_.begin(); v != other.variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        return mult;
    }

    AutoDf<ScalarType> operator-(AutoDf<ScalarType>& other)
    {
        AutoDf<ScalarType> mult(AutoType::kSubtractType);
        mult.first_ = this;
        mult.second_ = &other;
        *mult.value_ = *value_ - *other.value_;

        // copy all values together
        for (auto v = variables_.begin(); v != variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        // overwrites if some of them the same
        for (auto v = other.variables_.begin(); v != other.variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        return mult;
    }

    AutoDf<ScalarType> operator*(AutoDf<ScalarType>& other)
    {
        AutoDf<ScalarType> mult(AutoType::kMultType);
        mult.first_ = this;
        mult.second_ = &other;
        *mult.value_ = *value_ * *other.value_;

        // copy all values together
        for (auto v = variables_.begin(); v != variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        for (auto v = other.variables_.begin(); v != other.variables_.end(); v++)
        {
            mult.variables_[v->first] = v->second;
        }
        return mult;
    }

    AutoDf<ScalarType> operator/(AutoDf<ScalarType>& other)
    {
        AutoDf<ScalarType> divres(AutoType::kDivType);
        divres.first_ = this;
        divres.second_ = &other;
        *divres.value_ = *value_ / *other.value_;

        // copy all values together
        for (auto v = variables_.begin(); v != variables_.end(); v++)
        {
            divres.variables_[v->first] = v->second;
        }
        for (auto v = other.variables_.begin(); v != other.variables_.end(); v++)
        {
            divres.variables_[v->first] = v->second;
        }
        return divres;
    }

    ~AutoDf()
    {
        if (own_first)
        {
            delete first_;
            first_ = nullptr;
        }
        if (own_second)
        {
            delete second_;
            second_ = nullptr;
        }
        variables_.clear();
        value_.reset();
    }

    friend AutoDf<ScalarType> operator+(AutoDf<ScalarType>& scalar_value, const ScalarType scalar);
    friend AutoDf<ScalarType> operator+(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType> operator-(AutoDf<ScalarType>& scalar_value, const ScalarType scalar);
    friend AutoDf<ScalarType> operator-(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType> operator*(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType> operator*(AutoDf<ScalarType>& other, const ScalarType scalar_value);
    friend AutoDf<ScalarType> operator/(const ScalarType scalar_value, AutoDf<ScalarType>& other);
    friend AutoDf<ScalarType> operator/(AutoDf<ScalarType>& other, const ScalarType scalar_value);
};

AutoDf<float> operator+(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> sum(AutoDf<float>::AutoType::kSumType, nullptr, scalar);
    sum.first_ = &other;
    sum.variables_ = other.variables_;
    *sum.value_ = *other.value_ + scalar_value;
    return sum;
}

AutoDf<float> operator+(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> sum(AutoDf<float>::AutoType::kSumType, scalar, nullptr);
    sum.second_ = &other;
    sum.variables_ = other.variables_;
    *sum.value_ = *other.value_ + scalar_value;
    return sum;
}

AutoDf<float> operator-(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> sum(AutoDf<float>::AutoType::kScalarType, nullptr, scalar);
    sum.first_ = &other;
    sum.variables_ = other.variables_;
    *sum.value_ = *other.value_ - scalar_value;
    return sum;
}

AutoDf<float> operator-(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> sum(AutoDf<float>::AutoType::kScalarType, scalar, nullptr);
    sum.second_ = &other;
    sum.variables_ = other.variables_;
    *sum.value_ = *other.value_ - scalar_value;
    return sum;
}

AutoDf<float> operator*(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> mult(AutoDf<float>::AutoType::kMultType, scalar, nullptr);
    mult.second_ = &other;
    mult.variables_ = other.variables_;
    *mult.value_ = *other.value_ * scalar_value;
    return mult;
}

AutoDf<float> operator*(AutoDf<float>& other, const float scalar_value)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> mult(AutoDf<float>::AutoType::kMultType, nullptr, scalar);
    mult.first_ = &other;
    mult.variables_ = other.variables_;
    *mult.value_ = *other.value_ * scalar_value;
    return mult;
}

AutoDf<float> operator/(const float scalar_value, AutoDf<float>& other)
{
    AutoDf<float>* scalar = new AutoDf<float>(scalar_value, true);
    AutoDf<float> mult(AutoDf<float>::AutoType::kDivType, scalar, nullptr);
    mult.second_ = &other;
    mult.variables_ = other.variables_;
    *mult.value_ = scalar_value / *other.value_;
    return mult;
}

template <typename underlying = float>
AutoDf<float> operator/(AutoDf<underlying>& other, const underlying scalar_value)
{
    AutoDf<underlying>* scalar = new AutoDf<underlying>(scalar_value, true);
    AutoDf<underlying> mult(AutoDf<underlying>::AutoType::kDivType, nullptr, scalar);
    mult.first_ = &other;
    mult.variables_ = other.variables_;
    *mult.value_ = *other.value_ / scalar_value;
    return mult;
}

template <>
size_t AutoDf<float>::id_increment = 0U;

}  // namespace tiny_autodiff

#endif  // TINY_AUTODIFF_H_
