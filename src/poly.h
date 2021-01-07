// ----------------------------------------------------------------------------
//
// FILENAME: poly.h
//
// DESCRIPTION:
//    This file define the basic class, Poly, of this project. It contains some
//    basic operations of one polynomial.
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_POLY_H_
#define POLY_POLY_H_

#include <cassert>
#include <cmath>
#include <iostream>

#include "param.h"

// Represent a polynomial. The template "n" represent the maximum degree of this
// polynomial. Constant first.
//
// EXAMPLE:
//  x^2 + 2x + 1 can be in instantiate with class Poly<5>
//
// OPERATORS:
//  This class overload the some operators
//
//  -= : can accept anothor polynomial or a double
//  += : can accept anothor polynomial or a double
//  *= : Only accept double number
//  /= : Only accept double number
//  [] : Accept integer, return the corresponding coef
//  >> : Accept integer, right shift coefficients with n units
//
template <int n>
class Poly {
 public:
  // Constructor a zero polynomial
  inline Poly() : coef_{} { num_coef_ = 1; }

  // Construct a polynomial only idx-th degree has coefficient num
  inline Poly(const double num, const int idx) : coef_{} {
    coef_[idx] = num;
    num_coef_ = std::max(num_coef_, idx + 1);
  }

  // Construct a polynomial with specified coefficents.
  Poly(const double* input_coef, const int num_input) : coef_{} {
    assert(num_input <= n + 1);
    for (int i = 0; i < num_input; i++) coef_[i] = input_coef[i];
    num_coef_ = num_input;
  }

  // Copy construct
  template <int n1>
  inline Poly(const Poly<n1>& copy_poly) : coef_{} {
    static_assert(n1 <= n);
    for (int i = 0; i < n1; i++) coef_[i] = copy_poly[i];
    num_coef_ = copy_poly.Size();
  }

  // Copy assignment
  template <int n1>
  inline Poly<n>& operator=(const Poly<n1>& copy_poly) {
    static_assert(n1 <= n);
    for (int i = 0; i < n1; i++) coef_[i] = copy_poly[i];
    for (int i = n1; i <= n; i++) coef_[i] = 0.0;
    num_coef_ = copy_poly.Size();
    return *this;
  }

  // Move construct
  template <int n1>
  inline Poly(Poly<n1>&& move_poly) : coef_{} {
    static_assert(n1 <= n);
    for (int i = 0; i < n1; i++) coef_[i] = move_poly[i];
    num_coef_ = move_poly.Size();
  }

  // Move assignment
  template <int n1>
  inline Poly<n>& operator=(Poly<n1>&& move_poly) {
    static_assert(n1 <= n);
    for (int i = 0; i < n1; i++) coef_[i] = move_poly[i];
    for (int i = n1; i <= n; i++) coef_[i] = 0.0;
    num_coef_ = move_poly.Size();
    return *this;
  }

  // Get or set number of coefficents
  inline const int Size() const { return num_coef_; }
  inline void set_num_coef(int num_coef_input) { num_coef_ = num_coef_input; }

  // Get value of polynomial at point x
  inline double ValueAt(double x) {
    double ans = coef_[num_coef_ - 1];
    for (int i = num_coef_ - 2; i >= 0; i--) ans = ans * x + coef_[i];
    return ans;
  }

  // Get gradient of polynomial at point x
  inline double GradientAt(double x) {
    double ans = coef_[num_coef_ - 1] * (num_coef_ - 1);
    for (int i = num_coef_ - 2; i >= 1; i--) ans = ans * x + coef_[i] * i;
    return ans;
  }

  // Get leading coefficent
  inline double lead_coef() const { return coef_[num_coef_ - 1]; }

  // --------------------------------------------------------------------------
  //
  // Overload operators
  //
  // --------------------------------------------------------------------------

  // Get i-th coefficient
  inline const double& operator[](int i) const { return coef_[i]; }
  inline double& operator[](int i) { return coef_[i]; }

  // Get sum of this and poly2. Need to guarantee n1 is no more than n.
  // Then it will reset the num_coef_
  template <int n1>
  inline Poly<n>& operator+=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i < poly2.Size(); i++) coef_[i] += poly2[i];

    if (num_coef_ == poly2.Size())
      for (int i = num_coef_ - 1; i > 0 && std::abs(coef_[i]) < kEPSILON; i--)
        num_coef_--;
    else
      num_coef_ = std::max(num_coef_, poly2.Size());

    return *this;
  }

  // Sum of polynomial with a number
  inline Poly<n>& operator+=(const double num) {
    coef_[0] += num;
    return *this;
  }

  // Get difference of this and poly2. Need to guarantee n1 is no more than n.
  // Then it will reset the num_coef_ if two polynomial have same num of coef
  template <int n1>
  inline Poly<n>& operator-=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i < poly2.Size(); i++) coef_[i] -= poly2[i];

    if (num_coef_ == poly2.Size())
      for (int i = num_coef_ - 1; i > 0 && std::abs(coef_[i]) < kEPSILON; i--)
        num_coef_--;
    else
      num_coef_ = std::max(num_coef_, poly2.Size());

    return *this;
  }

  // Difference between polynomial and a number
  inline Poly<n>& operator-=(const double num) {
    coef_[0] -= num;
    return *this;
  }

  // Production of a polynomial and a number
  inline Poly<n>& operator*=(const double num) {
    for (int i = 0; i < num_coef_; i++) coef_[i] *= num;
    return *this;
  }

  // Division of a polynomial and a number
  inline Poly<n>& operator/=(const double num) {
    for (int i = 0; i < num_coef_; i++) coef_[i] /= num;
    return *this;
  }

  // Right shift operator, move coef right
  inline Poly<n>& operator>>(const int move_num) {
    if (move_num == 0) return *this;
    for (int i = num_coef_ - 1; i >= 0; i--) {
      coef_[i + move_num] = coef_[i];
      coef_[i] = 0.0;
    }
    num_coef_ += move_num;
    return *this;
  }

 private:
  double coef_[n + 1];
  int num_coef_;
};

// Struct for division result.
// Store both quotient and remainder
template <int n>
struct DivsionRet {
  Poly<n> quotient;
  Poly<n> remainder;
};

// Sum of two polynomials
template <int n1, int n2>
Poly<std::max(n1, n2)> operator+(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) += poly2;
  } else {
    return Poly<n2>(poly2) += poly1;
  }
}

// Difference of two polynomials
template <int n1, int n2>
Poly<std::max(n1, n2)> operator-(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) -= poly2;
  } else {
    return Poly<n2>(poly1) -= poly2;
  }
}

// Production of two polynomials.
template <int n1, int n2>
Poly<n1 + n2> operator*(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  Poly<n1 + n2> ret;
  ret.Size() = n1 + n2 + 1;
  for (int i = 0; i < poly1.Size(); i++)
    for (int j = 0; j < poly2.Size(); j++) ret[i + j] += poly1[i] * poly2[j];
  for (int i = n1 + n2; i > 0 && std::abs(ret[i]) < kEPSILON; i--) ret.Size()--;
  return ret;
}

// Production of polynomial and a number
template <int n>
Poly<n> operator*(const Poly<n>& poly, const double num) {
  return Poly<n>(poly) *= num;
}

// Production of polynomial and a number
template <int n>
Poly<n> operator*(const double num, const Poly<n>& poly) {
  return Poly<n>(poly) *= num;
}

// Division of two polynomials
// Applied long division method, pls check :
// https://en.wikipedia.org/wiki/Polynomial_long_division
template <int n1, int n2>
DivsionRet<n1> operator/(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  Poly<n1> quotient;
  Poly<n1> remainder(poly1);
  int degree = poly2.Size() - 1, remainder_degree = remainder.Size() - 1;
  double lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    double division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;
    quotient[degree_idx] += division;
    quotient.set_num_coef(std::max(quotient.Size(), degree_idx + 1));

    Poly<n1> sub(poly2);
    sub >> degree_idx;
    sub *= division;
    remainder -= sub;
    remainder_degree = remainder.Size() - 1;
  }
  return DivsionRet<n1>{quotient, remainder};
}

// Division of polynomial and a number
template <int n>
Poly<n> operator/(const Poly<n>& poly, const double num) {
  return Poly<n>(poly) /= num;
}

// Division of polynomial and a number
template <int n>
Poly<n> operator/(const double num, const Poly<n>& poly) {
  return Poly<n>(poly) /= num;
}

// Print out the polynomial. Ouput format will be like "cx^0 + bx^1 + ax^2"
//
// EXAMPLE:
//    Poly<3> a(Poly<3>(2, 3.0));
//    std::cout<<a<<std::endl;
//    It will print "3.0 x^3 "
template <int n>
static std::ostream& operator<<(std::ostream& out, const Poly<n>& u) {
  for (int i = 0; i < u.Size(); i++) {
    out << " ";
    if (u[i] >= 0) out << '+';
    out << u[i] << "x^" << i;
  }
  return out;
}

#endif  // POLY_POLYH_
