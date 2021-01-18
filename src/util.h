// ----------------------------------------------------------------------------
//
// FILENAME: util.h
//
// DESCRIPTION:
//    This file define some public funcitons that will be used in this project
//    These functions can be both used in Continued fraction and Budan method
//
// FUNCTIONS:
//    GCD                 : GCD of two polynomials
//    SquareFreeDecompose : Decompose polynomial into square-free
//    UpperBound          : Upper bound of roots
//    IsZero              : Return true if polynomial is zero
//    AddToX              : Change x to x+h
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

// TODO : const referenct
#ifndef POLY_UTIL_H
#define POLY_UTIL_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "poly.h"

// Return true if polynomial has no non-zero coefficient
template <int n>
bool IsZero(const Poly<n>& poly) {
  return poly.get_degree() == 0 && std::fabs(poly[0]) <= kEPSILON;
}

// This will return the Greatest Common Divider(GCD) of two polynomials.
// Since we do not need the remainder in this project, so it just return the
// GCD. Applied Euclid's Algorithm. ref
// https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor
//
// EXAMPLE:
//    Polynomial a = x^2 - 2x + 1
//    Polynomial b = x - 1
//    GCD(a, b) should return x - 1
//
//    // TODO : Too much copy
//    pass a poly to store result.
template <int n1, int n2>
Poly<std::min(n1, n2)> GCD(const Poly<n1>& poly1, const Poly<n2>& poly2) {
  if constexpr (n2 > n1)
    return GCD(poly2, poly1);
  else {
    auto divans = DivRemainder(poly1, poly2);
    if (IsZero(divans.remainder)) return poly2;
    Poly<n2> ret(GCD(poly2, divans.remainder));
    return ret;
  }
}

// This function will decompose a polynomial in to an array of square free
// polynomials. Applied Yun's algorithm, ref:
// https://en.wikipedia.org/wiki/Square-free_polynomial
// // TODO : avoid vector. using array of fixed polynomial. take array as
// arguments and return how many polys decompose out
template <int n>
std::vector<Poly<n>> SquareFreeDecompose(Poly<n>& poly) {
  std::vector<Poly<n>> ans;

  auto fd(poly.Derivative());
  auto a(GCD(poly, fd));
  auto b(DivRemainder(poly, a).quotient);
  auto c(DivRemainder(fd, a).quotient);
  auto d(c - b.Derivative());

  while (!(b.get_degree() == 0 && std::fabs(b[0] - 1) <= kEPSILON)) {  // b != 1
    if (IsZero(d)) {
      ans.emplace_back(b);
      break;
    }
    a = GCD(b, d);
    ans.emplace_back(a);
    b = DivRemainder(b, a).quotient;
    c = DivRemainder(d, a).quotient;
    d = c - b.Derivative();
  }

  return ans;
}

// Return the upper bound of value of roots
// Applied Cauchy's bound.
// p(x) = a0 + a1x + a2 x^2 ... + an x^n
// upperbound = 1 + max(abs(a0/an), abs(a1/an).... abs(an-1/an))
// https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
template <int n>
double UpperBound(const Poly<n>& poly) {
  double lc = poly.lead_coef(), ans = poly[0] / lc;
  for (int i = 1; i < poly.get_degree(); i++)
    ans = std::fmax(ans, std::fabs(poly[i] / lc));
  return 1 + ans;
}

// Change x to x + h. Applied Taylor Expansion to this,
// p(x+h) = p(h) + p'(h)x + 1/2*p''(h)x^2 ... 1/(n!) * p^n(h)*x^n
// time complexity: n^2
// ref: https://math.stackexchange.com/questions/694565/polynomial-shift
template <int n>
Poly<n> AddToX(const Poly<n>& poly, double h) {
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

#endif  // POLY_UTIL_H
