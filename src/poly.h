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
//  /= : Only accept double
//  *= : Only accept double
//  [] : Accept integer, return the corresponding coef
//
template <int n>
class Poly {
  static_assert(n >= 0);

 private:
  double coef_[n + 1];
  int degree_;

 public:
  // Constructor a zero polynomial
  Poly() : coef_{} { degree_ = 0; }

  // Construct a polynomial with specified coefficents.
  Poly(const double* input_coef, int num_input) : coef_{} {
    assert(num_input <= n + 1);
    for (int i = 0; i < num_input; i++) coef_[i] = input_coef[i];
    degree_ = num_input - 1;
  }

  // Copy from smaller Poly
  template <int n1>
  Poly(const Poly<n1>& copy_poly) : coef_{} {
    static_assert(n1 <= n);
    for (int i = 0; i <= n1; i++) coef_[i] = copy_poly[i];
    degree_ = copy_poly.get_degree();
  }

  // Assigned from smaller Poly
  template <int n1>
  Poly<n>& operator=(const Poly<n1>& copy_poly) {
    static_assert(n1 <= n);
    for (int i = 0; i <= n1; i++) coef_[i] = copy_poly[i];
    for (int i = n1 + 1; i <= n; i++) coef_[i] = 0.0;
    degree_ = copy_poly.get_degree();
    return *this;
  }

  const int get_degree() const { return degree_; }

  // Set degree, since it is defined by the coef_, so can not
  // directly changed from outside
  void set_degree() {
    degree_ = n;
    for (degree_ = n; degree_ > 0 && std::fabs(coef_[degree_]) < kEPSILON;
         degree_--)
      coef_[degree_] = 0.0;
  }

  // Get value of polynomial at point x
  const double ValueAt(double x) const {
    double ans = coef_[degree_];
    for (int i = degree_ - 1; i >= 0; i--) ans = ans * x + coef_[i];
    return ans;
  }

  // Get derivative of polynomial at point x
  const double DerivativeAt(double x) const {
    double ans = coef_[degree_] * (degree_);
    for (int i = degree_ - 1; i >= 1; i--) ans = ans * x + coef_[i] * i;
    return ans;
  }

  // Get derivative of polynomal
  const Poly<std::max(n - 1, 0)> Derivative() const {
    Poly<std::max(n - 1, 0)> ans;
    for (int i = 0; i < degree_; i++) ans[i] = coef_[i + 1] * (i + 1);
    ans.set_degree();
    return ans;
  }

  // Get leading coefficent
  const double lead_coef() const { return coef_[degree_]; }

  // Get sign change of coefficients
  const int SignChange() const {
    int ret = 0;
    bool prev = coef_[degree_] > 0;
    for (int i = degree_ - 1; i >= 0; i--) {
      if (std::fabs(coef_[i]) >= kEPSILON &&  // larger than 0
          (coef_[i] > 0 != prev)) {           // sign changed
        prev = !prev;
        ret++;
      }
    }
    return ret;
  }

  // --------------------------------------------------------------------------
  //
  // Overload operators
  //
  // --------------------------------------------------------------------------

  // Get or set i-th coefficient
  const double operator[](int i) const { return coef_[i]; }
  double& operator[](int i) { return coef_[i]; }

  // Get sum of this and poly2. Need to guarantee n1 is no more than n.
  // Then it will reset the degree
  template <int n1>
  const Poly<n>& operator+=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++) coef_[i] += poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  // Sum of polynomial with a number
  const Poly<n>& operator+=(double num) {
    coef_[0] += num;
    return *this;
  }

  // Get difference of this and poly2. Need to guarantee n1 is no more than n.
  // Then it will reset the num_coef_ if two polynomial have degree
  template <int n1>
  const Poly<n>& operator-=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++) coef_[i] -= poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  // Difference between polynomial and a number
  const Poly<n>& operator-=(double num) {
    coef_[0] -= num;
    return *this;
  }

  // Division of a polynomial and a number
  const Poly<n>& operator/=(double num) {
    for (int i = 0; i <= degree_; i++) coef_[i] /= num;
    set_degree();
    return *this;
  }

  // Production of a polynomial and a number
  const Poly<n>& operator*=(double num) {
    for (int i = 0; i < degree_; i++) coef_[i] *= num;
    set_degree();
    return *this;
  }
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
  for (int i = 0; i <= poly1.get_degree(); i++)
    for (int j = 0; j <= poly1.get_degree(); j++)
      ret[i + j] += poly1[i] * poly2[j];
  ret.set_degree();
  return ret;
}

// Production of polynomial and a number
template <int n>
Poly<n> operator*(const Poly<n>& poly, double num) {
  return Poly<n>(poly) *= num;
}

// Production of polynomial and a number
template <int n>
Poly<n> operator*(double num, const Poly<n>& poly) {
  return Poly<n>(poly) *= num;
}

// Struct for DivRemainder result.
// Store both quotient and remainder
template <int n1, int n2>
struct DivsionRet {
  Poly<n1> quotient;
  Poly<n2> remainder;
};

// poly1 -= (poly2 >> move_num) * scale
// poly1 will be modified in this function
// Since this only used in division, it is guaranteed that n2 + move_num <= n1
template <int n1, int n2>
void MinusRightMoveScale(const Poly<n2>& poly2, int move_num, double scale,
                         Poly<n1>& poly1) {
  for (int i = poly2.get_degree(); i >= 0; i--) {
    poly1[i + move_num] -= poly2[i] * scale;
  }
  poly1.set_degree();
}

// Division of two polynomials
// Applied long division method, pls check :
// https://en.wikipedia.org/wiki/Polynomial_long_division
template <int n1, int n2>
DivsionRet<n1, std::max(n2 - 1, 0)> DivRemainder(const Poly<n1>& poly1,
                                                 const Poly<n2>& poly2) {
  DivsionRet<n1, std::max(n2 - 1, 0)> ret;

  // If poly2 is a constant number
  if (poly2.get_degree() == 0) {
    ret.quotient = poly1 / poly2[0];
    return ret;
  }

  Poly<n1> remainder(poly1);
  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  double lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    double division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;
    ret.quotient[degree_idx] += division;

    MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }
  ret.quotient.set_degree();

  // Set remaineder in return struct
  for (int i = 0; i <= remainder.get_degree(); i++)
    ret.remainder[i] = remainder[i];
  ret.remainder.set_degree();

  return ret;
}

// Division of polynomial and a number
template <int n>
Poly<n> operator/(const Poly<n>& poly, double num) {
  return Poly<n>(poly) /= num;
}

//
// EXAMPLE:
//    Poly<3> a(Poly<3>(2, 3.0));
//    std::cout<<a<<std::endl;
//    It will print "3.0 x^3 "
template <int n>
static std::ostream& operator<<(std::ostream& out, const Poly<n>& u) {
  for (int i = 0; i <= u.get_degree(); i++) {
    out << " ";
    if (u[i] >= 0) out << '+';
    out << u[i] << "x^" << i;
  }
  return out;
}

#endif  // POLY_POLYH_
