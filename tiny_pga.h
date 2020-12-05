//
// C++ Implementation of Projectove Geometric Algebra (a.k.a Plane-based Geometric Algebra)
// Highly-templatized implementation helps with number of issues:
// * reducing computational complexity of PGA operators via compile-time optimizations
// * optimizing memory footprint
// * checking correctness of multivector assignments / blades matching
//
// Created by Sergey Smirnov on 5.12.2020.
// Email: sergei.smirnov@gmail.com
// Copyright Sergey Smirnov  / Seregium Oy 2020
//

#ifndef TINY_PGA_H
#define TINY_PGA_H

#include <array>
#include <cstdint>
#include <type_traits>

namespace tiny_pga
{
using Elems = std::uint16_t;

namespace elems
{

enum class BitValues : Elems
{
  kScalar = (1U << 0),
  kE0 = (1U << 1),
  kE1 = (1U << 2),
  kE2 = (1U << 3),
  kE3 = (1U << 4),
  kE01 = (1U << 5),
  kE02 = (1U << 6),
  kE03 = (1U << 7),
  kE12 = (1U << 8),
  kE31 = (1U << 9),
  kE23 = (1U << 10),
  kE021 = (1U << 11),
  kE013 = (1U << 12),
  kE032 = (1U << 13),
  kE123 = (1U << 14),
  kE0123 = (1U << 15)
};

constexpr bool has_scalar(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kScalar);
}

constexpr bool has_e0(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE0);
}

constexpr bool has_e1(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE1);
}

constexpr bool has_e2(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE2);
}

constexpr bool has_e3(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE3);
}

constexpr bool has_e01(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE01);
}

constexpr bool has_e02(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE02);
}

constexpr bool has_e03(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE03);
}

constexpr bool has_e12(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE12);
}

constexpr bool has_e31(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE31);
}

constexpr bool has_e23(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE23);
}

constexpr bool has_e021(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE021);
}

constexpr bool has_e013(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE013);
}

constexpr bool has_e032(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE032);
}

constexpr bool has_e123(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE123);
}

constexpr bool has_e0123(const Elems elem)
{
  return elem & static_cast<Elems>(BitValues::kE0123);
  ;
}

constexpr bool has_vector(Elems elems)
{
  return has_e0(elems) || has_e1(elems) || has_e2(elems) || has_e3(elems);
}

constexpr bool has_bivectorE(Elems elems)
{
  return has_e23(elems) || has_e31(elems) || has_e12(elems) || has_scalar(elems);
}

constexpr bool has_bivector0(Elems elems)
{
  return has_e01(elems) || has_e02(elems) || has_e03(elems) || has_e0123(elems);
}

constexpr bool has_trivector(Elems elems)
{
  return has_e021(elems) || has_e013(elems) || has_e032(elems) || has_e123(elems);
}

constexpr Elems multiplication(const Elems elems1, const Elems elems2)
{
  Elems out_elements = 0U;

  const bool scalar =
      (has_scalar(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e1(elems2)) ||
      (has_e2(elems1) && has_e2(elems2)) || (has_e3(elems1) && has_e3(elems2)) ||
      (has_e12(elems1) && has_e12(elems2)) || (has_e31(elems1) && has_e31(elems2)) ||
      (has_e23(elems1) && has_e23(elems2)) || (has_e123(elems1) && has_e123(elems2));

  const bool e0 =
      (has_scalar(elems1) && has_e0(elems2)) || (has_e0(elems1) && has_scalar(elems2)) ||
      (has_e1(elems1) && has_e01(elems2)) || (has_e2(elems1) && has_e02(elems2)) ||
      (has_e3(elems1) && has_e02(elems2)) || (has_e01(elems1) && has_e1(elems2)) ||
      (has_e02(elems1) && has_e2(elems2)) || (has_e03(elems1) && has_e3(elems2)) ||
      (has_e12(elems1) && has_e021(elems2)) || (has_e31(elems1) && has_e013(elems2)) ||
      (has_e23(elems1) && has_e032(elems2)) || (has_e021(elems1) && has_e12(elems2)) ||
      (has_e013(elems1) && has_e31(elems2)) || (has_e032(elems1) && has_e23(elems2)) ||
      (has_e123(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_e123(elems2));

  const bool e1 = (has_scalar(elems1) && has_e1(elems2)) ||
                  (has_e1(elems1) && has_scalar(elems2)) || (has_e2(elems1) && has_e12(elems2)) ||
                  (has_e3(elems1) && has_e31(elems2)) || (has_e12(elems1) && has_e2(elems2)) ||
                  (has_e31(elems1) && has_e3(elems2)) || (has_e23(elems1) && has_e123(elems2)) ||
                  (has_e123(elems1) && has_e23(elems2));

  const bool e2 = (has_scalar(elems1) && has_e2(elems2)) ||
                  (has_e2(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e12(elems2)) ||
                  (has_e3(elems1) && has_e23(elems2)) || (has_e12(elems1) && has_e1(elems2)) ||
                  (has_e23(elems1) && has_e3(elems2)) || (has_e31(elems1) && has_e123(elems2)) ||
                  (has_e123(elems1) && has_e31(elems2));

  const bool e3 = (has_scalar(elems1) && has_e3(elems2)) ||
                  (has_e3(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e31(elems2)) ||
                  (has_e2(elems1) && has_e23(elems2)) || (has_e31(elems1) && has_e1(elems2)) ||
                  (has_e23(elems1) && has_e2(elems2)) || (has_e12(elems1) && has_e123(elems2)) ||
                  (has_e123(elems1) && has_e12(elems2));

  const bool e01 = (has_scalar(elems1) && has_e01(elems2)) ||
                   (has_e01(elems1) && has_scalar(elems2)) || (has_e0(elems1) && has_e1(elems2)) ||
                   (has_e1(elems1) && has_e0(elems2)) || (has_e2(elems1) && has_e021(elems2)) ||
                   (has_e3(elems1) && has_e013(elems2)) || (has_e021(elems1) && has_e2(elems2)) ||
                   (has_e013(elems1) && has_e3(elems2)) || (has_e02(elems1) && has_e12(elems2)) ||
                   (has_e03(elems1) && has_e31(elems2)) || (has_e12(elems1) && has_e02(elems2)) ||
                   (has_e31(elems1) && has_e03(elems2)) || (has_e032(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e032(elems2));

  const bool e02 = (has_scalar(elems1) && has_e02(elems2)) ||
                   (has_e02(elems1) && has_scalar(elems2)) || (has_e0(elems1) && has_e2(elems2)) ||
                   (has_e2(elems1) && has_e0(elems2)) || (has_e1(elems1) && has_e021(elems2)) ||
                   (has_e3(elems1) && has_e032(elems2)) || (has_e021(elems1) && has_e1(elems2)) ||
                   (has_e032(elems1) && has_e3(elems2)) || (has_e01(elems1) && has_e12(elems2)) ||
                   (has_e03(elems1) && has_e23(elems2)) || (has_e12(elems1) && has_e01(elems2)) ||
                   (has_e23(elems1) && has_e03(elems2)) || (has_e013(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e013(elems2));

  const bool e03 = (has_scalar(elems1) && has_e03(elems2)) ||
                   (has_e03(elems1) && has_scalar(elems2)) || (has_e0(elems1) && has_e3(elems2)) ||
                   (has_e3(elems1) && has_e0(elems2)) || (has_e1(elems1) && has_e013(elems2)) ||
                   (has_e2(elems1) && has_e032(elems2)) || (has_e013(elems1) && has_e1(elems2)) ||
                   (has_e032(elems1) && has_e2(elems2)) || (has_e01(elems1) && has_e31(elems2)) ||
                   (has_e02(elems1) && has_e23(elems2)) || (has_e31(elems1) && has_e01(elems2)) ||
                   (has_e23(elems1) && has_e02(elems2)) || (has_e021(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e021(elems2));

  const bool e12 = (has_scalar(elems1) && has_e12(elems2)) ||
                   (has_e12(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e2(elems2)) ||
                   (has_e2(elems1) && has_e1(elems2)) || (has_e3(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e3(elems2)) || (has_e31(elems1) && has_e23(elems2)) ||
                   (has_e23(elems1) && has_e31(elems2));

  const bool e31 = (has_scalar(elems1) && has_e31(elems2)) ||
                   (has_e31(elems1) && has_scalar(elems2)) || (has_e1(elems1) && has_e3(elems2)) ||
                   (has_e3(elems1) && has_e1(elems2)) || (has_e2(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e2(elems2)) || (has_e12(elems1) && has_e23(elems2)) ||
                   (has_e23(elems1) && has_e12(elems2));

  const bool e23 = (has_scalar(elems1) && has_e23(elems2)) ||
                   (has_e23(elems1) && has_scalar(elems2)) || (has_e2(elems1) && has_e3(elems2)) ||
                   (has_e3(elems1) && has_e2(elems2)) || (has_e1(elems1) && has_e123(elems2)) ||
                   (has_e123(elems1) && has_e1(elems2)) || (has_e12(elems1) && has_e31(elems2)) ||
                   (has_e31(elems1) && has_e12(elems2));

  const bool e021 =
      (has_scalar(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_scalar(elems2)) ||
      (has_e0(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e0(elems2)) ||
      (has_e1(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e1(elems2)) ||
      (has_e2(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e2(elems2)) ||
      (has_e3(elems1) && has_e0123(elems2)) || (has_e123(elems1) && has_e3(elems2)) ||
      (has_e03(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e03(elems2)) ||
      (has_e23(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e23(elems2));

  const bool e013 =
      (has_scalar(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_scalar(elems2)) ||
      (has_e0(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e0(elems2)) ||
      (has_e1(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e1(elems2)) ||
      (has_e3(elems1) && has_e01(elems2)) || (has_e01(elems1) && has_e3(elems2)) ||
      (has_e2(elems1) && has_e0123(elems2)) || (has_e123(elems1) && has_e2(elems2)) ||
      (has_e02(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e02(elems2)) ||
      (has_e23(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e23(elems2));

  const bool e032 =
      (has_scalar(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_scalar(elems2)) ||
      (has_e0(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e0(elems2)) ||
      (has_e2(elems1) && has_e03(elems2)) || (has_e03(elems1) && has_e2(elems2)) ||
      (has_e3(elems1) && has_e02(elems2)) || (has_e02(elems1) && has_e3(elems2)) ||
      (has_e1(elems1) && has_e0123(elems2)) || (has_e123(elems1) && has_e1(elems2)) ||
      (has_e01(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e01(elems2)) ||
      (has_e31(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e31(elems2));

  const bool e123 = (has_scalar(elems1) && has_e123(elems2)) ||
                    (has_e123(elems1) && has_scalar(elems2)) ||
                    (has_e1(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e1(elems2)) ||
                    (has_e2(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e2(elems2)) ||
                    (has_e3(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e3(elems2));

  const bool e0123 =
      (has_scalar(elems1) && has_e0123(elems2)) || (has_e0123(elems1) && has_scalar(elems2)) ||
      (has_e0(elems1) && has_e123(elems2)) || (has_e123(elems1) && has_e0(elems2)) ||
      (has_e1(elems1) && has_e032(elems2)) || (has_e032(elems1) && has_e1(elems2)) ||
      (has_e2(elems1) && has_e013(elems2)) || (has_e013(elems1) && has_e02(elems2)) ||
      (has_e3(elems1) && has_e021(elems2)) || (has_e021(elems1) && has_e3(elems2)) ||
      (has_e01(elems1) && has_e23(elems2)) || (has_e23(elems1) && has_e01(elems2)) ||
      (has_e02(elems1) && has_e31(elems2)) || (has_e31(elems1) && has_e02(elems2)) ||
      (has_e03(elems1) && has_e12(elems2)) || (has_e12(elems1) && has_e03(elems2));

  out_elements |= scalar ? static_cast<Elems>(BitValues::kScalar) : 0;
  out_elements |= e0 ? static_cast<Elems>(BitValues::kE0) : 0;
  out_elements |= e1 ? static_cast<Elems>(BitValues::kE1) : 0;
  out_elements |= e2 ? static_cast<Elems>(BitValues::kE2) : 0;
  out_elements |= e3 ? static_cast<Elems>(BitValues::kE3) : 0;
  out_elements |= e01 ? static_cast<Elems>(BitValues::kE01) : 0;
  out_elements |= e02 ? static_cast<Elems>(BitValues::kE02) : 0;
  out_elements |= e03 ? static_cast<Elems>(BitValues::kE03) : 0;
  out_elements |= e12 ? static_cast<Elems>(BitValues::kE12) : 0;
  out_elements |= e31 ? static_cast<Elems>(BitValues::kE31) : 0;
  out_elements |= e23 ? static_cast<Elems>(BitValues::kE23) : 0;
  out_elements |= e021 ? static_cast<Elems>(BitValues::kE021) : 0;
  out_elements |= e013 ? static_cast<Elems>(BitValues::kE013) : 0;
  out_elements |= e032 ? static_cast<Elems>(BitValues::kE032) : 0;
  out_elements |= e123 ? static_cast<Elems>(BitValues::kE123) : 0;
  out_elements |= e0123 ? static_cast<Elems>(BitValues::kE0123) : 0;

  return out_elements;
}

constexpr Elems PlaneElems =
    static_cast<Elems>(elems::BitValues::kE0) | static_cast<Elems>(elems::BitValues::kE1) |
    static_cast<Elems>(elems::BitValues::kE2) | static_cast<Elems>(elems::BitValues::kE3);

constexpr Elems LineElems =
    static_cast<Elems>(elems::BitValues::kE01) | static_cast<Elems>(elems::BitValues::kE02) |
    static_cast<Elems>(elems::BitValues::kE03) | static_cast<Elems>(elems::BitValues::kE23) |
    static_cast<Elems>(elems::BitValues::kE31) | static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems PointElems =
    static_cast<Elems>(elems::BitValues::kE123) | static_cast<Elems>(elems::BitValues::kE032) |
    static_cast<Elems>(elems::BitValues::kE013) | static_cast<Elems>(elems::BitValues::kE021);

constexpr Elems RotorElems =
    static_cast<Elems>(elems::BitValues::kScalar) | static_cast<Elems>(elems::BitValues::kE23) |
    static_cast<Elems>(elems::BitValues::kE31) | static_cast<Elems>(elems::BitValues::kE12);

constexpr Elems TranslatorElems = static_cast<Elems>(elems::BitValues::kE01) |
                                  static_cast<Elems>(elems::BitValues::kE02) |
                                  static_cast<Elems>(elems::BitValues::kE03);

constexpr Elems MotorElems =
    static_cast<Elems>(elems::BitValues::kScalar) | static_cast<Elems>(elems::BitValues::kE23) |
    static_cast<Elems>(elems::BitValues::kE31) | static_cast<Elems>(elems::BitValues::kE12) |
    static_cast<Elems>(elems::BitValues::kE01) | static_cast<Elems>(elems::BitValues::kE02) |
    static_cast<Elems>(elems::BitValues::kE03) | static_cast<Elems>(elems::BitValues::kE0123);

} // namespace elems

template <Elems elements> struct Multivector
{
  template <bool Condition, typename T> struct Conditional
  {
    T value;
  };
  template <typename T> struct Conditional<false, T>
  {
  };

  // Optimization of Memory Footprint with use of conditional elements
  Conditional<elems::has_vector(elements), std::array<float, 4U>> Vector;
  Conditional<elems::has_bivectorE(elements), std::array<float, 4U>> BivectorE;
  Conditional<elems::has_bivector0(elements), std::array<float, 4U>> Bivector0;
  Conditional<elems::has_trivector(elements), std::array<float, 4U>> Trivector;

  template <typename T = float>
  typename std::enable_if<elems::has_scalar(elements), T>::type &scalar()
  {
    return BivectorE.value[0];
  }

  /// Alias to @ref e0123() function
  template <typename T = float>
  typename std::enable_if<elems::has_e0123(elements), T>::type &pseudo()
  {
    return Bivector0.value[3];
  }

  template <typename T = float> typename std::enable_if<elems::has_e0(elements), T>::type &e0()
  {
    return Vector.value[0];
  }

  template <typename T = float> typename std::enable_if<elems::has_e1(elements), T>::type &e1()
  {
    return Vector.value[1];
  }

  template <typename T = float> typename std::enable_if<elems::has_e2(elements), T>::type &e2()
  {
    return Vector.value[2];
  }

  template <typename T = float> typename std::enable_if<elems::has_e3(elements), T>::type &e3()
  {
    return Vector.value[3];
  }

  template <typename T = float> typename std::enable_if<elems::has_e01(elements), T>::type &e01()
  {
    return Bivector0.value[0];
  }

  template <typename T = float> typename std::enable_if<elems::has_e02(elements), T>::type &e02()
  {
    return Bivector0.value[1];
  }
  template <typename T = float> typename std::enable_if<elems::has_e03(elements), T>::type &e03()
  {
    return Bivector0.value[2];
  }

  ///
  template <typename T = float>
  typename std::enable_if<elems::has_e0123(elements), T>::type &e0123()
  {
    return Bivector0.value[3];
  }

  template <typename T = float> typename std::enable_if<elems::has_e12(elements), T>::type &e12()
  {
    return BivectorE.value[1];
  }

  template <typename T = float> typename std::enable_if<elems::has_e31(elements), T>::type &e31()
  {
    return BivectorE.value[2];
  }

  template <typename T = float> typename std::enable_if<elems::has_e23(elements), T>::type &e23()
  {
    return BivectorE.value[3];
  }

  template <typename T = float> typename std::enable_if<elems::has_e021(elements), T>::type &e021()
  {
    return Trivector.value[0];
  }

  template <typename T = float> typename std::enable_if<elems::has_e013(elements), T>::type &e013()
  {
    return Trivector.value[1];
  }

  template <typename T = float> typename std::enable_if<elems::has_e032(elements), T>::type &e032()
  {
    return Trivector.value[2];
  }

  template <typename T = float> typename std::enable_if<elems::has_e123(elements), T>::type &e123()
  {
    return Trivector.value[2];
  }

  template <class T = Multivector<elems::RotorElems>>
  typename std::enable_if<elems::has_bivectorE(elements), T>::type rotor()
  {
    Multivector<elems::RotorElems> rotor{BivectorE.value};
    return rotor;
  }

  template <class T = Multivector<elems::TranslatorElems>>
  typename std::enable_if<elems::has_bivector0(elements), T>::type translator()
  {
    Multivector<elems::TranslatorElems> translator{Bivector0.value};
    return translator;
  }

  template <class T = Multivector<elems::MotorElems>>
  typename std::enable_if<elems::has_bivector0(elements), T>::type motor()
  {
    Multivector<elems::MotorElems> motor{BivectorE.value, Bivector0.value};
    return motor;
  }

  template <Elems other_elements>
  Multivector<elems::multiplication(elements, other_elements)> operator*(
      const Multivector<other_elements> &other)
  {
    constexpr Elems out_elems = elems::multiplication(elements, other_elements);
    Multivector<out_elems> out{};

    if (elems::has_scalar(out_elems))
    {
      if (elems::has_scalar(elements) && elems::has_scalar(other_elements))
      {
        out.scalar() += scalar() * other.scalar();
      }

      if (elems::has_e1(elements) && elems::has_e1(other_elements))
      {
        out.scalar() += e1() * other.e1();
      }

      if (elems::has_e2(elements) && elems::has_e2(other_elements))
      {
        out.scalar() += e2() * other.e2();
      }

      if (elems::has_e3(elements) && elems::has_e3(other_elements))
      {
        out.scalar() += e3() * other.e3();
      }

      if (elems::has_e12(elements) && elems::has_e12(other_elements))
      {
        out.scalar() -= e12() * other.e12();
      }

      if (elems::has_e31(elements) && elems::has_e31(other_elements))
      {
        out.scalar() -= e31() * other.e31();
      }

      if (elems::has_e23(elements) && elems::has_e23(other_elements))
      {
        out.scalar() -= e23() * other.e23();
      }

      if (elems::has_e123(elements) && elems::has_e123(other_elements))
      {
        out.scalar() -= e123() * other.e123();
      }
    }

    if (elems::has_e0(out_elems))
    {
      if (elems::has_scalar(elements) && elems::has_e0(other_elements))
      {
        out.e0() += scalar() * other.e0();
      }

      if (elems::has_e0(elements) && elems::has_scalar(other_elements))
      {
        out.e0() += e0() * other.scalar();
      }

      if (elems::has_e1(elements) && elems::has_e01(other_elements))
      {
        out.e0() -= e1() * other.e01();
      }

      if (elems::has_e01(elements) && elems::has_e1(other_elements))
      {
        out.e0() += e01() * other.e1();
      }

      if (elems::has_e2(elements) && elems::has_e02(other_elements))
      {
        out.e0() -= e2() * other.e02();
      }

      if (elems::has_e02(elements) && elems::has_e2(other_elements))
      {
        out.e0() += e02() * other.e2();
      }

      if (elems::has_e3(elements) && elems::has_e03(other_elements))
      {
        out.e0() -= e3() * other.e03();
      }

      if (elems::has_e03(elements) && elems::has_e3(other_elements))
      {
        out.e0() += e03() * other.e3();
      }

      if (elems::has_e021(elements) && elems::has_e12(other_elements))
      {
        out.e0() += e021() * other.e12();
      }

      if (elems::has_e12(elements) && elems::has_e021(other_elements))
      {
        out.e0() += e12() * other.e021();
      }

      if (elems::has_e31(elements) && elems::has_e013(other_elements))
      {
        out.e0() += e31() * other.e013();
      }

      if (elems::has_e23(elements) && elems::has_e032(other_elements))
      {
        out.e0() += e23() * other.e032();
      }

      if (elems::has_e013(elements) && elems::has_e31(other_elements))
      {
        out.e0() += e013() * other.e31();
      }

      if (elems::has_e032(elements) && elems::has_e23(other_elements))
      {
        out.e0() += e032() * other.e23();
      }

      if (elems::has_e123(elements) && elems::has_e0123(other_elements))
      {
        out.e0() += e123() * other.pseudo();
      }

      if (elems::has_e0123(elements) && elems::has_e123(other_elements))
      {
        out.e0() += pseudo() * other.e123();
      }
    }

    return out;
  }
};

using Plane = Multivector<elems::PlaneElems>;
using Line = Multivector<elems::LineElems>;
using Point = Multivector<elems::PointElems>;
using Rotor = Multivector<elems::RotorElems>;
using Translator = Multivector<elems::TranslatorElems>;
using Motor = Multivector<elems::MotorElems>;

} // namespace tiny_pga

#endif // TINY_PGA_H
