// ----------------------------------------------------------------------------
//
// FILENAME: util.h
//
// DESCRIPTION:
//    This file define some public funcitons that will be used in this project
//    These functions can be both used in Continued fraction and Budan method
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_UTIL_H
#define POLY_UTIL_H

#include <algorithm>
#include <cmath>

#include "poly.h"
#include "range.h"

/**
 * @brief :return True if all coefficient of poly is zero
 */
template <int n> bool IsZero(const Poly<n> &poly) {
  return poly.get_degree() == 0 && std::fabs(poly[0]) <= kEPSILON;
}

/**
 * @brief Calculate GCD of poly1 and poly2, and save it into ret
 * Applied Euclid's Algorithm. ref:
 * https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
 *
 * @tparam n1 :maximum degree of poly1
 * @tparam n2 :maximum degree of poly2
 * @tparam n3 :maximum degree of ret
 * @param poly1 :polynomial 1
 * @param poly2 :polynomial 2
 * @param ret :store GCD of poly1 and poly2
 */
template <int n1, int n2, int n3>
void GCD_helper_(const Poly<n1> &poly1, const Poly<n2> &poly2, Poly<n3> &ret) {
  static_assert(n3 >= n2);
  auto remainder = Remainder(poly1, poly2);
  if (IsZero(remainder)) {
    ret = poly2;
    return;
  }
  GCD_helper_(poly2, remainder, ret);
}

/**
 * @brief :Warp GCD_helper_ function return the GCD of two polynomials
 *
 * @tparam n1 :maximum degree of poly1
 * @tparam n2 :maximum degree of poly2
 * @param poly1 :polynomial 1
 * @param poly2 :polynomial 2
 */
template <int n1, int n2>
Poly<std::min(n1, n2)> GCD(const Poly<n1> &poly1, const Poly<n2> &poly2) {
  if constexpr (n2 > n1)
    return GCD(poly2, poly1);
  else {
    Poly<n2> ret;
    GCD_helper_(poly1, poly2, ret);
    return ret;
  }
}

/**
 * @brief Decopose a polynomial into an array of square free polynomials
 * Applied Yun's algrithm, ref:
 * https://en.wikipedia.org/wiki/Square-free_polynomial
 *
 * @tparam n :maximum degree of polynomal
 * @param poly : polynomial
 * @param ans :array of polynomials used to restore the result
 * @return : number of square free polynomials
 */
template <int n> int SquareFreeDecompose(const Poly<n> &poly, Poly<n> *ans) {
  int ret = 0; // number of square free polynomial

  auto fd(poly.Derivative());
  auto a(GCD(poly, fd));
  auto b(Quotient(poly, a));
  auto c(Quotient(fd, a));
  auto d(c - b.Derivative());

  while (!(b.get_degree() == 0 && std::fabs(b[0] - 1) <= kEPSILON)) { // b != 1
    if (IsZero(d)) {
      ans[ret++] = b;
      break;
    }
    a = GCD(b, d);
    ans[ret++] = a;
    b = Quotient(b, a);
    c = Quotient(d, a);
    d = c - b.Derivative();
  }
  return ret;
}

/**
 * @brief Calculate the upper bound of value of roots
 * Applied Cauchy's bound.
 * p(x) = a0 + a1x + a2 x^2 ... + an x^n
 * upperbound = 1 + max(abs(a0/an), abs(a1/an).... abs(an-1/an))
 * ref:
 * https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
 *
 * @tparam n :maximum degree of polynomial
 * @param poly :polynomial
 * @return :upper bound of roots
 */
template <int n> double UpperBound(const Poly<n> &poly) {
  double lc = poly.lead_coef(), ans = std::fabs(poly[0] / lc);
  for (int i = 1; i < poly.get_degree(); i++)
    ans = std::fmax(ans, std::fabs(poly[i] / lc));
  return 1 + ans;
}

/**
 * @brief Replace "x" in poly with "x+h"
 * Applied Taylor Expansion to this.
 * p(x+h) = p(h) + p'(h)x + 1/2*p''(h)x^2 ... 1/(n!) * p^n(h)*x^n
 * ref: https://math.stackexchange.com/questions/694565/polynomial-shift
 *
 * @tparam n :maximum degree of poly
 * @param poly :polynomial
 * @param h :number add to x
 * @return :polynomial that replace x in poly to x+h
 */
template <int n> Poly<n> AddToX(const Poly<n> &poly, double h) {
  if (h == 0)
    return poly;
  Poly<n> ret, tmp(poly);
  ret[0] = tmp.ValueAt(h);
  double divisor = 1;

  for (int i = 1; i <= poly.get_degree(); i++) {
    tmp = tmp.Derivative();
    divisor *= i;
    ret[i] = tmp.ValueAt(h) * 1 / divisor;
  }

  ret.set_degree();
  return ret;
}

/**
 * @brief :Check if zero is root of polynomial. If it is, remove it. Since It is
 * square free polynomial. So it at most have one zero root.
 *
 * @tparam n :Maximum degree of poly
 * @param repeat_time :Repeat time of this polynomial in original polynomial
 * @param poly :Polynomial, might be modified.
 * @param ranges :Store results, might be modified
 * @param num_roots :Sotore the number of roots, might be modified
 */
template <int n>
void HandleZeroRoots(int repeat_time, Poly<n> *poly, Range *ranges,
                     int *num_roots) {
  // Zero is not root
  if (std::fabs((*poly)[0]) >= kEPSILON)
    return;

  // remove zero root
  for (int i = 1; i <= poly->get_degree(); i++)
    (*poly)[i - 1] = (*poly)[i];
  (*poly)[poly->get_degree()] = 0.0;
  poly->set_degree(poly->get_degree() - 1);

  // add to ranges
  Range ans{0.0, 0.0};
  for (int i = 0; i < repeat_time; i++) {
    ranges[*num_roots] = ans;
    (*num_roots)++;
  }

  return;
}

/**
 * @brief :Hanle linear polynomial, It has no zero root
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param repeat_time :Repeat time of this polynomial in original polynomial
 * @param ranges :Store results, might be modified
 * @param num_roots :Sotore the number of roots, might be modified
 */
template <int n>
void HandleLinear(const Poly<n> &poly, int repeat_time, Range *ranges,
                  int *num_roots) {
  double root(-poly[0] / poly[1]);
  Range ans{root, root};
  for (int i = 0; i < repeat_time; i++) {
    ranges[*num_roots] = ans;
    (*num_roots)++;
  }

  return;
}

/**
 * @brief :Hanle quadratic polynomial, It has no zero root
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param repeat_time :Repeat time of this polynomial in original polynomial
 * @param ranges :Store results, might be modified
 * @param num_roots :Sotore the number of roots, might be modified
 */
template <int n>
void HandleQuadratic(const Poly<n> &poly, int repeat_time, Range *ranges,
                     int *num_roots) {
  double delta(poly[1] * poly[1] - 4 * poly[0] * poly[2]); // b^2 - 4ac

  if (delta <= 0)
    return;
  delta = std::sqrt(delta);

  Range ans1, ans2;
  ans1.left_end = ans1.right_end = (-poly[1] + delta) / (2 * poly[2]);
  ans2.left_end = ans2.right_end = (-poly[1] - delta) / (2 * poly[2]);

  for (int i = 0; i < repeat_time; i++) {
    ranges[*num_roots] = ans1;
    ranges[(*num_roots) + repeat_time] = ans2;
    (*num_roots)++;
  }
  (*num_roots) += repeat_time;

  return;
}

#endif // POLY_UTIL_H
