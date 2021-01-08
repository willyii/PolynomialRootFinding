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
//    GCD
//
// TODO:
//    SquareFreeDecompose
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

  return Poly<n1>(GCD(poly2, divans.remainder));
}

#endif  // POLY_UTIL_H
