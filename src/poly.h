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

/**
 * @brief This class used to represent a polynomial with maximum degree n.
 *
 * @tparam n: maximum degree of this polynomial.
 */
template <int n>
class Poly {
  static_assert(n >= 0);

 private:
  double coef_[n + 1];
  int degree_;

 public:
  Poly() : coef_{} { degree_ = 0; }

  /**
   * @brief contruct polynomial with sepcific coefficients.
   *
   * @param input_coef: pointer to coefficients array.
   * @param num_input: number of coefficients
   */
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

  int get_degree() const { return degree_; }

  void set_degree() {
    for (degree_ = n; degree_ > 0 && std::fabs(coef_[degree_]) < kEPSILON;
         degree_--)
      coef_[degree_] = 0.0;
  }
  void set_degree(int degree_input) {
    assert(degree_input <= n);
    degree_ = degree_input;
  }

  // Get value of polynomial at point x
  double ValueAt(double x) const {
    double ans = coef_[degree_];
    for (int i = degree_ - 1; i >= 0; i--) ans = ans * x + coef_[i];
    return ans;
  }

  // Get derivative of polynomial at point x
  double DerivativeAt(double x) const {
    double ans = coef_[degree_] * degree_;
    for (int i = degree_ - 1; i >= 1; i--) ans = ans * x + coef_[i] * i;
    return ans;
  }

  // Get derivative of polynomal
  Poly<std::max(n - 1, 0)> Derivative() const {
    Poly<std::max(n - 1, 0)> ans;
    for (int i = 0; i < degree_; i++) ans[i] = coef_[i + 1] * (i + 1);
    ans.set_degree(degree_ - 1);
    return ans;
  }

  // Get leading coefficent
  double lead_coef() const { return coef_[degree_]; }

  // Get # of sign change of coefficients
  int SignChange() const {
    int ret = 0;
    bool prev = coef_[degree_] > 0;
    for (int i = degree_ - 1; i >= 0; i--) {
      if (std::fabs(coef_[i]) >= kEPSILON &&  // not 0
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

  double operator[](int i) const { return coef_[i]; }
  double& operator[](int i) { return coef_[i]; }

  template <int n1>
  Poly<n>& operator+=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++) coef_[i] += poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  Poly<n>& operator+=(double num) {
    coef_[0] += num;
    return *this;
  }

  template <int n1>
  Poly<n>& operator-=(const Poly<n1>& poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++) coef_[i] -= poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  Poly<n>& operator-=(double num) {
    coef_[0] -= num;
    return *this;
  }

  Poly<n>& operator/=(double num) {
    for (int i = 0; i <= degree_; i++) coef_[i] /= num;
    return *this;
  }

  Poly<n>& operator*=(double num) {
    for (int i = 0; i < degree_; i++) coef_[i] *= num;
    return *this;
  }
};

template <int n1, int n2>
Poly<std::max(n1, n2)> operator+(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) += poly2;
  } else {
    return Poly<n2>(poly2) += poly1;
  }
}

template <int n1, int n2>
Poly<std::max(n1, n2)> operator-(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) -= poly2;
  } else {
    return Poly<n2>(poly1) -= poly2;
  }
}

template <int n1, int n2>
Poly<n1 + n2> operator*(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  Poly<n1 + n2> ret;
  for (int i = 0; i <= poly1.get_degree(); i++)
    for (int j = 0; j <= poly1.get_degree(); j++)
      ret[i + j] += poly1[i] * poly2[j];
  ret.set_degree(poly1.get_degree() - poly2.get_degree());
  return ret;
}

template <int n>
Poly<n> operator*(const Poly<n>& poly, double num) {
  return Poly<n>(poly) *= num;
}

template <int n>
Poly<n> operator*(double num, const Poly<n>& poly) {
  return Poly<n>(poly) *= num;
}

/**
 * struct DivsionRet - division result
 *
 * @tparam n1: maximum degree of quotient
 * @tparam n2: maximum degree of remainder
 */
template <int n1, int n2>
struct DivsionRet {
  Poly<n1> quotient;
  Poly<n2> remainder;
};

/**
 * @brief This function will calculate poly1 -= (poly2 >> move_num) * scale.
 * This function only used in division, It's guaranteed that poly2.degree +
 * num_num <= n1
 *
 * @tparam n1: maximum degree of poly1
 * @tparam n2: maximum degree of poly2
 * @param poly2: polynomial 2
 * @param move_num: number of unit that coef of poly2 will moved
 * @param scale: scale coefficent of poly2 with scale
 * @param poly1: polynomial 1
 */
template <int n1, int n2>
void MinusRightMoveScale(const Poly<n2>& poly2, int move_num, double scale,
                         Poly<n1>& poly1) {
  for (int i = poly2.get_degree(); i >= 0; i--) {
    poly1[i + move_num] -= poly2[i] * scale;
  }
  poly1.set_degree();
}

/**
 * @brief poly1/poly2. Return both quotient and remainder
 * Applied Euclidean division, ref :
 * https://en.wikipedia.org/wiki/Polynomial_long_division
 *
 * @tparam n1: maximum degree of polynomial 1
 * @tparam n2: maximum degree of polynomial 2
 * @param poly1: polynomial 1
 * @param poly2: polynomial 2
 */
template <int n1, int n2>
DivsionRet<n1, std::max(n2 - 1, 0)> Division(const Poly<n1>& poly1,
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
    ret.quotient[degree_idx] = division;

    MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }
  ret.quotient.set_degree(poly1.get_degree() - poly2.get_degree());

  // Set remaineder in return struct
  for (int i = 0; i <= remainder.get_degree(); i++)
    ret.remainder[i] = remainder[i];
  ret.remainder.set_degree(remainder.get_degree());

  return ret;
}

/**
 * @brief Quotient of poly1/poly2
 *
 * @tparam n1 :maximum degree of poly1
 * @tparam n2 :maximum degree of poly2
 * @param poly1 :polynomial 1
 * @param poly2 :polynomial 2
 */
template <int n1, int n2>
Poly<n1> Quotient(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  Poly<n1> quotient, remainder(poly1);

  // If poly2 is a constant number
  if (poly2.get_degree() == 0) {
    quotient = poly1 / poly2[0];
    return quotient;
  }

  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  double lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    double division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;
    quotient[degree_idx] = division;

    MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }
  quotient.set_degree(poly1.get_degree() - poly2.get_degree());

  return quotient;
}

/**
 * @brief :Remainder of poly1/poly2
 *
 * @tparam n1 :maximum degree of poly1
 * @tparam n2 :maximum degree of poly2
 * @param poly1 :polynomial 1
 * @param poly2 :polynomial 2
 */
template <int n1, int n2>
Poly<std::max(n2 - 1, 0)> Remainder(const Poly<n1>& poly1,
                                    const Poly<n2>& poly2) {
  Poly<n1> remainder(poly1);
  Poly<std::max(n2 - 1, 0)> remainder_ret;

  // If poly2 is a constant number
  if (poly2.get_degree() == 0) return remainder_ret;

  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  double lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    double division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;

    MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }

  // Set remaineder should returned
  for (int i = 0; i <= remainder.get_degree(); i++)
    remainder_ret[i] = remainder[i];
  remainder_ret.set_degree(remainder.get_degree());

  return remainder_ret;
}

template <int n>
Poly<n> operator/(const Poly<n>& poly, double num) {
  return Poly<n>(poly) /= num;
}

// Print out polynomial u
template <int n>
static std::ostream& operator<<(std::ostream& out, const Poly<n>& u) {
  for (int i = 0; i <= u.get_degree(); i++) {
    out << " ";
    if (std::fabs(u[i]) < kEPSILON) continue;
    if (u[i] >= 0) out << '+';
    out << u[i] << "x^" << i;
  }
  return out;
}

#endif  // POLY_POLYH_
