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
//    Replace             : Replace the x with polynomial
//    NumSignChange       : Number of sign change among coefficients
//    UpperBound          : Upper bound of roots
//    LowerBound          : Lower bound of roots
//    IsZero              : Return true if polynomial is zero
//
// DOING:
//
// TODO:
//    GCD # Test
//    SquareFreeDecompose # Test
//    Replace
//    NumSignChange
//    UpperBound
//    LowerBound
//
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_UTIL_H
#define POLY_UTIL_H

#include <vector>

#include "poly.h"

// Return true if polynomial has no non-zero coefficient
template <int n>
bool IsZero(Poly<n> poly) {
  return poly.Size() == 1 && std::fabs(poly[0]) <= kEPSILON;
}

// This will return the Greatest Common Divider(GCD) of two polynomials.
// Since we do not need the remainder in this project, so it just return the
// GCD.
//
// EXAMPLE:
//    Polynomial a = x^2 - 2x + 1
//    Polynomial b = x - 1
//    GCD(a, b) should return x - 1
template <int n1, int n2>
Poly<n1> GCD(Poly<n1>& poly1, Poly<n2>& poly2) {
  if (IsZero(poly2)) return poly1;
  auto divans = poly1 / poly2;
  return Poly<n1>(std::move(GCD(poly2, divans.remainder)));
}

// This function will decompose a polynomial in to an array of square free
// polynomials.
template <int n>
std::vector<Poly<n>> SquareFreeDecompose(Poly<n> poly) {
  std::vector<Poly<n>> ans;

  auto fd(poly.Derivative());
  auto a(GCD(poly, fd));
  auto b((poly / a).quotient);
  auto c((fd / a).quotient);
  auto d(c - b.Derivative());

  while (!(b.Size() == 1 && std::fabs(b[0] - 1) <= kEPSILON)) {  // b != 1
    a = GCD(b, d);
    ans.emplace_back(a);
    b = (b / a).quotient;
    c = (d / a).quotient;
    d = c - b.Derivative();
  }

  return ans;
}

#endif  // POLY_UTIL_H
