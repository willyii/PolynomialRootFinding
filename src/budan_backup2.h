// ----------------------------------------------------------------------------
//
// FILENAME: budan.h
//
// DESCRIPTION:
//    This file define and implement functions that will be used in root
//    isolation with budan theorem.
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_BUDAN_H
#define POLY_BUDAN_H

#include <cmath>

#include "poly.h"
#include "range.h"
#include "util.h"

// debug
#include <stdio.h>

/**
 * Isolate the root of a square free polynomial, and store it to ranges.
 * Applied bisection method
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param duplicate_times :Repeat time of polynomial in original polynomial
 * @param left :Left end of search range
 * @param left_change : Number of sign change of p(x+left_end)
 * @param right :Right end of search range
 * @param right_change : Number of sign change of p(x+right_end)
 * @param ranges :Store isolation results, might be modified
 * @param num_roots : Store the number of roots, might be modified
 */
template <int n>
void BudanSquareFreeSolve(const Poly<n> &poly, int duplicate_times, double left,
                          int left_change, double right, int right_change,
                          Range *ranges, int *num_roots) {
  // no root in this range
  if (left_change == right_change)
    return;

  // search range smaller than threshold
  else if ((right - left) < kMINRANGE) {
    // else if ((right - left) < kMINRANGE || (left_change - right_change) == 1)
    // {
    AddToRange((left_change - right_change), left, right, ranges, num_roots);
    return;
  }

  double mid((left + right) * 0.5);
  // if (std::fabs(poly.ValueAt(mid)) < kEPSILON) {
  //  AddToRange(duplicate_times, mid, mid, ranges, num_roots);
  //  return;
  //}
  int mid_change(AddToX(poly, mid).SignChange());

  // left side
  BudanSquareFreeSolve<n>(poly, duplicate_times, left, left_change, mid,
                          mid_change, ranges, num_roots);
  // right side
  BudanSquareFreeSolve<n>(poly, duplicate_times, mid, mid_change, right,
                          right_change, ranges, num_roots);
}

/**
 * Apply Budan' Theorem to get range of every root of polynomial that
 * represented by coef and coef_num. Duplicated roots will represented by
 * same range.
 *
 * @param coef :Pointer to coefficients of polynomial
 * @param coef_num :Number of coefficients
 * @param ranges :Store isolation results, might be modified
 * @return :Number of roots
 */
int BudanRootIsolate(const double *coef, int coef_num, Range *ranges) {

  Poly<kMAXDEGREE> original_poly(coef, coef_num);

  std::cout << "DEBUG: original poly " << original_poly << std::endl;

  int num_roots(0);

  // Handle zero roots
  AddToRange(ZeroRoots(&original_poly), 0.0, 0.0, ranges, &num_roots);

  if (original_poly.get_degree() == 0) // constant
    return num_roots;
  else if (original_poly.get_degree() == 1) // linear
    Linear<kMAXDEGREE>(original_poly, 1, ranges, &num_roots);
  else if (original_poly.get_degree() == 2) // quadratic
    Quadratic<kMAXDEGREE>(original_poly, 1, ranges, &num_roots);
  else {
    double right(UpperBound(original_poly)), left = -right;
    int left_change(AddToX(original_poly, left).SignChange());
    int right_change(AddToX(original_poly, right).SignChange());

    BudanSquareFreeSolve<kMAXDEGREE>(original_poly, 1, left, left_change, right,
                                     right_change, ranges, &num_roots);
  }

  return num_roots;
}

#endif // POLY_BUDAN_H
