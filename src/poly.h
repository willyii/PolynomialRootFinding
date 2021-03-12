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

#include <algorithm>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/utility_fwd.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h>

#include "param.h"

typedef boost::numeric::interval<double> interval;

/**
 * Class represent polynomials with maximum degree n. Actual degree might be
 * smaller than n
 *
 * @tparam n :Maximum degree of this polynomial.
 */
template <int n> class Poly {
  static_assert(n >= 0);

private:
  interval coef_[n + 1];
  int degree_;

public:
  Poly() : coef_{} { degree_ = 0; }

  /**
   * Contruct polynomial with sepcific coefficients.
   *
   * @param input_coef: Pointer to coefficients array.
   * @param num_input: Number of coefficients
   */
  Poly(const double *input_coef, int num_input) : coef_{} {

    assert(num_input <= n + 1);
    for (int i = 0; i < num_input; i++) {
      coef_[i] = input_coef[i];
    }
    set_degree();
  }

  /**
   * Copy from smaller Polynomial
   */
  template <int n1> Poly(const Poly<n1> &copy_poly) : coef_{} {
    static_assert(n1 <= n);
    for (int i = 0; i <= n1; i++)
      coef_[i] = copy_poly[i];
    degree_ = copy_poly.get_degree();
  }

  /**
   * Assign from smaller Polynomial
   */
  template <int n1> Poly<n> &operator=(const Poly<n1> &copy_poly) {
    static_assert(n1 <= n);
    for (int i = 0; i <= n1; i++)
      coef_[i] = copy_poly[i];
    for (int i = n1 + 1; i <= n; i++)
      coef_[i] = 0.0;
    degree_ = copy_poly.get_degree();
    return *this;
  }

  int get_degree() const { return degree_; }

  void set_degree() {
    for (degree_ = n; degree_ > 0 && containZero(degree_); degree_--)
      ;
  }

  void set_degree(int degree_input) {
    assert(degree_input <= n);
    degree_ = degree_input;
  }

  /**
   * Value of polynomial at point x
   */
  interval ValueAt(interval x) const {
    interval ans = coef_[degree_];
    for (int i = degree_ - 1; i >= 0; i--)
      ans = ans * x + coef_[i];
    return ans;
  }

  /**
   * Derivative of polynomial at point x
   */
  interval DerivativeAt(interval x) const {
    interval ans = coef_[degree_] * double(degree_);
    for (int i = degree_ - 1; i >= 1; i--)
      ans = ans * x + coef_[i] * double(i);
    return ans;
  }

  /**
   * Derivative of polynomal
   */
  Poly<std::max(n - 1, 0)> Derivative() const {
    Poly<std::max(n - 1, 0)> ans;
    for (int i = 0; i < degree_; i++)
      ans[i] = coef_[i + 1] * double(i + 1);
    for (int i = degree_; i <= n; i++) {
      ans[i] = 0;
    }
    ans.set_degree(degree_ - 1);
    return ans;
  }

  /**
   * Derivative of polynomal
   */
  void Derivative_() {
    for (int i = 0; i < degree_; i++)
      coef_[i] = coef_[i + 1] * double(i + 1);
    for (int i = degree_; i <= n; i++) {
      coef_[i] = 0;
    }
    degree_--;
  }

  interval lead_coef() const { return coef_[degree_]; }

  /**
   * Get number of sign changes in coefficients
   */
  int SignChange() const {
    int ret = 0;
    bool prev = coef_[degree_].lower() > 0;
    for (int i = degree_ - 1; i >= 0; i--) {
      if (!containZero(i) && (coef_[i].lower() > 0 != prev)) { // sign changed
        prev = !prev;
        ret++;
      }
    }
    return ret;
  }

  bool containZero(int i) const { return boost::numeric::zero_in(coef_[i]); }

  /**
   * --------------------------------------------------------------------------
   *
   * Overload operators
   *
   * --------------------------------------------------------------------------
   */

  interval operator[](int i) const { return coef_[i]; }
  interval &operator[](int i) { return coef_[i]; }

  template <int n1> Poly<n> &operator+=(const Poly<n1> &poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++)
      coef_[i] += poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  Poly<n> &operator+=(interval num) {
    coef_[0] += num;
    return *this;
  }

  template <int n1> Poly<n> &operator-=(const Poly<n1> &poly2) {
    static_assert(n1 <= n);
    for (int i = 0; i <= poly2.get_degree(); i++)
      coef_[i] -= poly2[i];

    if (degree_ == poly2.get_degree())
      set_degree();
    else
      degree_ = std::max(degree_, poly2.get_degree());

    return *this;
  }

  Poly<n> &operator-=(interval num) {
    coef_[0] -= num;
    return *this;
  }

  Poly<n> &operator*=(interval num) {
    for (int i = 0; i <= degree_; i++)
      coef_[i] *= num;
    return *this;
  }

  Poly<n> &operator/=(interval num) {
    for (int i = 0; i <= degree_; i++)
      coef_[i] /= num;
    return *this;
  }
};

/*
 * + operator
 */
template <int n1, int n2>
Poly<std::max(n1, n2)> operator+(const Poly<n1> &poly1, const Poly<n2> &poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) += poly2;
  } else {
    return Poly<n2>(poly2) += poly1;
  }
}

template <int n> Poly<n> operator+(const Poly<n> &poly, interval num) {
  return Poly<n>(poly) += num;
}

template <int n> Poly<n> operator+(interval num, const Poly<n> &poly) {
  return Poly<n>(poly) += num;
}

/*
 * - operator
 */
template <int n1, int n2>
Poly<std::max(n1, n2)> operator-(const Poly<n1> &poly1, const Poly<n2> &poly2) {
  if constexpr (n1 >= n2) {
    return Poly<n1>(poly1) -= poly2;
  } else {
    return Poly<n2>(poly1) -= poly2;
  }
}

template <int n> Poly<n> operator-(const Poly<n> &poly, interval num) {
  return Poly<n>(poly) -= num;
}

template <int n> Poly<n> operator-(interval num, const Poly<n> &poly) {
  return Poly<n>(poly * -1.0) += num;
}

/*
 * * operator
 */
template <int n1, int n2>
Poly<n1 + n2> operator*(const Poly<n1> &poly1, const Poly<n2> &poly2) {
  Poly<n1 + n2> ret;
  for (int i = 0; i <= poly1.get_degree(); i++)
    for (int j = 0; j <= poly1.get_degree(); j++)
      ret[i + j] += poly1[i] * poly2[j];
  ret.set_degree(poly1.get_degree() + poly2.get_degree());
  return ret;
}

template <int n> Poly<n> operator*(const Poly<n> &poly, interval num) {
  return Poly<n>(poly) *= num;
}

template <int n> Poly<n> operator*(interval num, const Poly<n> &poly) {
  return Poly<n>(poly) *= num;
}

/**
 * struct DivsionRet - division result
 *
 * @tparam n1: maximum degree of quotient
 * @tparam n2: maximum degree of remainder
 */
template <int n1, int n2> struct DivsionRet {
  Poly<n1> quotient;  /* Division quotient */
  Poly<n2> remainder; /* Division remainder */
};

/**
 * Helper funtion in divison. Calculate poly1 -= (poly2 >> move_num) * scale.
 * Since its only used in division, It's guaranteed that poly2.degree +
 * num_num <= n1
 *
 * @tparam n1 :Maximum degree of poly1
 * @tparam n2 :Maximum degree of poly2
 * @param poly2 :Polynomial 2
 * @param move_num :Number of unit poly2 will be moved
 * @param scale :Scale coefficent of poly2
 * @param poly1 :Polynomial 1, might be modified
 */
template <int n1, int n2>
void MinusRightMoveScale(const Poly<n2> &poly2, int move_num, interval scale,
                         Poly<n1> &poly1) {
  for (int i = poly2.get_degree(); i >= 0; i--) {
    poly1[i + move_num] -= poly2[i] * scale / poly2.lead_coef();
  }
  poly1.set_degree();
}

/**
 * poly1/poly2. Return both quotient and remainder
 * Applied Euclidean division, ref :
 * https://en.wikipedia.org/wiki/Polynomial_long_division
 *
 * @tparam n1 :Maximum degree of polynomial 1
 * @tparam n2 :Maximum degree of polynomial 2
 * @param poly1 :Polynomial 1
 * @param poly2 :Polynomial 2
 */
template <int n1, int n2>
DivsionRet<n1, std::max(n2 - 1, 0)> Division(const Poly<n1> &poly1,
                                             const Poly<n2> &poly2) {
  DivsionRet<n1, std::max(n2 - 1, 0)> ret;

  // If poly2 is a constant number
  if (poly2.get_degree() == 0) {
    ret.quotient = poly1 / poly2[0];
    return ret;
  }

  Poly<n1> remainder(poly1);
  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  interval lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    interval division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;
    ret.quotient[degree_idx] = division;

    MinusRightMoveScale(poly2, degree_idx, remainder.lead_coef(), remainder);
    // MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }
  ret.quotient.set_degree();

  // Set remaineder in return struct
  for (int i = 0; i <= remainder.get_degree(); i++)
    ret.remainder[i] = remainder[i];
  ret.remainder.set_degree(remainder.get_degree());

  return ret;
}

/**
 * Quotient of poly1/poly2
 *
 * @tparam n1 :Maximum degree of poly1
 * @tparam n2 :Maximum degree of poly2
 * @param poly1 :Polynomial 1
 * @param poly2 :Polynomial 2
 */
template <int n1, int n2>
Poly<n1> Quotient(const Poly<n1> &poly1, const Poly<n2> &poly2) {
  // std::cout << poly1 << " \ndiv\n " << poly2 << std::endl;
  Poly<n1> quotient, remainder(poly1);

  // If poly2 is a constant number
  if (poly2.get_degree() == 0) {
    quotient = poly1 / poly2[0];
    return quotient;
  }

  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  interval lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    interval division = remainder.lead_coef() / lead_coef;
    int degree_idx = remainder_degree - degree;
    quotient[degree_idx] = division;

    MinusRightMoveScale(poly2, degree_idx, remainder.lead_coef(), remainder);
    // MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();
  }

  quotient.set_degree();
  return quotient;
}

/**
 * Remainder of poly1/poly2
 *
 * @tparam n1 :Maximum degree of poly1
 * @tparam n2 :Maximum degree of poly2
 * @param poly1 :Polynomial 1
 * @param poly2 :Polynomial 2
 */
template <int n1, int n2>
Poly<std::max(n2 - 1, 0)> Remainder(const Poly<n1> &poly1,
                                    const Poly<n2> &poly2) {

  Poly<n1> remainder(poly1);
  Poly<std::max(n2 - 1, 0)> remainder_ret;

  // If poly2 is a constant number
  if (poly2.get_degree() == 0)
    return remainder_ret;

  int degree = poly2.get_degree(), remainder_degree = remainder.get_degree();
  interval lead_coef = poly2.lead_coef();

  while (remainder_degree >= degree) {
    interval division = remainder.lead_coef() / lead_coef;

    // std::cout << "Debug Remainder:  division "
    //<< boost::numeric::median(division) << "["
    //<< boost::numeric::width(division) << "]" << std::endl;
    int degree_idx = remainder_degree - degree;
    MinusRightMoveScale(poly2, degree_idx, remainder.lead_coef(), remainder);
    // MinusRightMoveScale(poly2, degree_idx, division, remainder);
    remainder_degree = remainder.get_degree();

    // std::cout << "Debug Remainder:  remainder " << remainder << std::endl;
  }

  // Set remaineder should returned
  for (int i = 0; i <= remainder.get_degree(); i++)
    remainder_ret[i] = remainder[i];
  remainder_ret.set_degree(remainder.get_degree());

  return remainder_ret;
}

template <int n> Poly<n> operator/(const Poly<n> &poly, interval num) {
  return Poly<n>(poly) /= num;
}

/**
 * Print out polynomial u
 */
template <int n>
static std::ostream &operator<<(std::ostream &out, const Poly<n> &u) {
  for (int i = 0; i <= u.get_degree(); i++) {
    // if (u.containZero(i))
    //  continue;
    if (boost::numeric::median(u[i]) >= 0)
      out << '+';
    // out << boost::numeric::median(u[i]) << "*x^" << i;
    out << boost::numeric::median(u[i]) << "[" << u[i].upper() - u[i].lower()
        << "]*x^" << i;
  }
  out << " degree = " << u.get_degree();
  return out;
}

#endif // POLY_POLYH_
