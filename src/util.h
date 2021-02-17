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
#include <boost/numeric/interval/utility_fwd.hpp>
#include <cassert>
#include <cmath>
#include <math.h>
#include <tuple>

#include "param.h"
#include "poly.h"
#include "range.h"

/**
 * Return True if all coefficients of poly is zero
 */
template <int n> bool IsZero(const Poly<n> &poly) {
  for (int i = poly.get_degree(); i >= 0; i--) {
    if (!poly.containZero(i)) {
      return false;
    }
  }
  return true;
}

template <int n> bool EndGCD(const Poly<n> &poly, int spacial) {
  for (int i = 0; i <= poly.get_degree(); i++) {
    double mid = boost::numeric::median(poly[i]);
    double tolerance = (poly[i].upper() - poly[i].lower()) / 2.0;
    tolerance *= std::pow(10, spacial);
    if (!(mid - tolerance <= 0.0 && mid + tolerance >= 0.0) &&
        std::abs(mid) > 1e-8)
      return false;
  }
  return true;
  return IsZero(poly);
}

/**
 * Make a polynomial monic
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial, might be modified.
 */
template <int n> void Monic(Poly<n> &poly) {

  // interval max_coef = poly[0];
  // for (int i = 1; i <= poly.get_degree(); i++)
  // max_coef = boost::numeric::max(max_coef, poly[i]);
  interval max_coef = poly.lead_coef();
  for (int i = 0; i <= poly.get_degree(); i++)
    poly[i] /= max_coef;
}

/**
 * Calculate GCD of poly1 and poly2, and save it into ret
 * Applied Euclid's Algorithm. ref:
 * https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
 *
 * @tparam n1 :Maximum degree of poly1
 * @tparam n2 :Maximum degree of poly2
 * @tparam n3 :Maximum degree of ret
 * @param poly1 :Polynomial 1
 * @param poly2 :Polynomial 2
 * @param ret :Store GCD of poly1 and poly2, might be modified
 */
template <int n1, int n2, int n3>
void GCD_helper_(const Poly<n1> &poly1, const Poly<n2> &poly2, Poly<n3> &ret,
                 int spacial) {
  static_assert(n3 >= n2);
  auto remainder = Remainder(poly1, poly2);
  std::cout << "- DEGBUG_GCD: remainder " << remainder << std::endl;
  if (EndGCD(remainder, spacial)) {
    ret = poly2;
    return;
  }
  GCD_helper_(poly2, remainder, ret, spacial);
}

/**
 * Warp GCD_helper_ function return the GCD of two polynomials
 *
 * @tparam n1 :Maximum degree of poly1
 * @tparam n2 :Maximum degree of poly2
 * @param poly1 :Polynomial 1
 * @param poly2 :Polynomial 2
 * @return :GCD of two polynomials
 */
template <int n1, int n2>
Poly<std::min(n1, n2)> GCD(const Poly<n1> &poly1, const Poly<n2> &poly2,
                           int spacial) {
  std::cout << "\n========GCD===========\n"
            << poly1 << "\n"
            << poly2 << std::endl;
  if constexpr (n2 > n1)
    return GCD(poly2, poly1, spacial);
  else {
    Poly<n2> ret;
    spacial = tmp(poly1);
    GCD_helper_(poly1, poly2, ret, spacial);
    Monic(ret);
    std::cout << "\n ========== END GCD ========= \n " << std::endl;
    return ret;
  }
}

template <int n> int tmp(const Poly<n> &poly) {
  double ret = std::abs(boost::numeric::median(poly[0]));
  for (int i = 1; i < poly.get_degree(); i++)
    ret = std::max(ret, std::abs(boost::numeric::median(poly[i])));
  return std::max(0, int(std::log10(std::abs(ret))) - 1);
}

/**
 * Decopose a polynomial into an array of square free polynomials
 * Applied Yun's algrithm, ref:
 * https://en.wikipedia.org/wiki/Square-free_polynomial
 *
 * @tparam n :Maximum degree of polynomal
 * @param poly :Polynomial to be decomposed
 * @param ans :Pointer to array of polynomials to restore the result, might be
 *             modifed
 * @return :Number of square free polynomials
 */
template <int n> int SquareFreeDecompose(const Poly<n> &poly_in, Poly<n> *ans) {
  Poly<n> poly(poly_in);
  Monic(poly);
  int ret = 0; // number of square free polynomial
  int end_degree = 0;

  int spacial = tmp(poly);
  std::cout << " FUCK : " << spacial << std::endl;
  auto fd(poly.Derivative());
  auto a(GCD(poly, fd, spacial));
  // if (a.get_degree() == 0) {
  //  ans[ret++] = poly;
  //  return ret;
  //}
  auto b(Quotient(poly, a));
  end_degree += b.get_degree();
  auto c(Quotient(fd, a));
  auto d(c - b.Derivative());

  std::cout << "DEBUG: a " << a << std::endl;

  std::cout << "DEBUG: b " << b << std::endl;
  std::cout << "DEBUG: c " << c << std::endl;
  std::cout << "DEBUG: d " << d << std::endl;
  std::cout << "================" << std::endl;
  while (!(b.get_degree() == 0)) { // b !=1
    if (end_degree >= poly.get_degree()) {
      ans[ret++] = b;
      assert(ret <= kMAXDEGREE);
      break;
    }
    // if (IsZero(d)) {
    // ans[ret++] = b;
    // assert(ret <= kMAXDEGREE);
    // break;
    //}
    a = GCD(b, d, spacial);

    ans[ret++] = a;
    assert(ret <= kMAXDEGREE);
    b = Quotient(b, a);
    end_degree += b.get_degree();
    c = Quotient(d, a);
    d = c - b.Derivative();
    std::cout << "DEBUG: a " << a << std::endl;
    std::cout << "DEBUG: b " << b << std::endl;
    std::cout << "DEBUG: c " << c << std::endl;
    std::cout << "DEBUG: d " << d << std::endl;
    std::cout << "================" << std::endl;
  }
  return ret;
}

/**
 * Calculate the upper bound of roots
 * Applied Cauchy's bound.
 * p(x) = a0 + a1x + a2 x^2 ... + an x^n
 * upperbound = 1 + max(abs(a0/an), abs(a1/an).... abs(an-1/an))
 * ref:
 * https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @return :Upper bound of roots
 */
template <int n> interval UpperBound(const Poly<n> &poly) {
  interval lc = poly.lead_coef(), ans = boost::numeric::abs(poly[0] / lc);
  for (int i = 1; i < poly.get_degree(); i++)
    ans = boost::numeric::max(ans, boost::numeric::abs(poly[i] / lc));
  return 1.0 + ans;
}

/**
 * Replace "x" in polynomial with "x+h"
 * Applied Taylor Expansion to this.
 * p(x+h) = p(h) + p'(h)x + 1/2*p''(h)x^2 ... 1/(n!) * p^n(h)*x^n
 * ref: https://math.stackexchange.com/questions/694565/polynomial-shift
 *
 * @tparam n :Maximum degree of poly
 * @param poly :Polynomial
 * @param h :Number add to x
 * @return :Polynomial that replace x in poly to x+h
 */
template <int n> Poly<n> AddToX(const Poly<n> &poly, interval h) {
  if (h.lower() <= 0.0 && h.upper() >= 0.0)
    return poly;
  Poly<n> ret, tmp(poly);
  ret[0] = tmp.ValueAt(h);
  interval divisor(1);

  for (int i = 1; i <= poly.get_degree(); i++) {
    tmp = tmp.Derivative();
    divisor *= i;
    ret[i] = tmp.ValueAt(h) / divisor;
  }

  ret.set_degree();
  return ret;
}

/**
 * Add range represented by left and right to ans
 *
 * @param repeat_time :Repeat time of this root
 * @param left :Left end of this root range
 * @param right :Right end of this root range
 * @param ranges :Store isolation results, might be modified
 * @param num_roots :Store the number of roots, might be modified
 */
/* TODO : median */
void AddToRange(int repeat_time, interval left, interval right, Range *ranges,
                int *num_roots) {
  if (repeat_time == 0)
    return;
  for (int i = 0; i < repeat_time; i++) {
    ranges[*num_roots].left_end = std::min(right.lower(), left.lower());
    ranges[*num_roots].right_end = std::max(right.upper(), left.upper());
    (*num_roots)++;
  }
  return;
}

/**
 * Return if zero is root of poly and remove zero root of polynomial. Since
 * its square free, there at most 1 zero root. This function won't add zero to
 * final ans;
 *
 * @tparam n :Maxiumu degree of polynomial
 * @param poly :Polynomial, might be modified
 * @return :True if zero is root of polynomial
 */
template <int n> bool ZeroRoots(Poly<n> *poly) {
  // Zero is not root
  if (!(*poly).containZero(0))
    return false;

  // remove zero root
  for (int i = 1; i <= poly->get_degree(); i++)
    (*poly)[i - 1] = (*poly)[i];
  (*poly)[poly->get_degree()] = 0.0;
  poly->set_degree(std::max(poly->get_degree() - 1, 0));

  return true;
}

/**
 * Solve linear polynomial, It has no zero root
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param repeat_time :Repeat time of this polynomial in original polynomial
 * @param ranges :Store results, might be modified
 * @param num_roots :Sotore the number of roots, might be modified
 */
template <int n>
void Linear(const Poly<n> &poly, int repeat_time, Range *ranges,
            int *num_roots) {
  interval root((-poly[0] / poly[1]));

  AddToRange(repeat_time, root, root, ranges, num_roots);
  return;
}

/**
 * Solve quadratic polynomial, It has no zero root
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param repeat_time :Repeat time of this polynomial in original polynomial
 * @param ranges :Store results, might be modified
 * @param num_roots :Sotore the number of roots, might be modified
 */
template <int n>
void Quadratic(const Poly<n> &poly, int repeat_time, Range *ranges,
               int *num_roots) {
  interval delta((poly[1] * poly[1] - 4.0 * poly[0] * poly[2])); // b^2 - 4ac

  if ((delta.upper() < 0))
    return;
  delta = boost::numeric::sqrt(delta);

  interval root1(-((poly[1] + delta) / (2.0 * poly[2])));
  AddToRange(repeat_time, root1, root1, ranges, num_roots);
  interval root2(-((poly[1] - delta) / (2.0 * poly[2])));
  AddToRange(repeat_time, root2, root2, ranges, num_roots);
  return;
}

#endif // POLY_UTIL_H
