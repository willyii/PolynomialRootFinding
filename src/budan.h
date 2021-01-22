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

/**
 * @brief :Isolate the root of a square free polynomial, Applied bisection
 * method
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param duplicate_times :How many times thie poly repeat in original poly
 * @param left : Left end of search range
 * @param left_change : Number of sign change in left end
 * @param right : Right end of search range
 * @param right_change : Number of sign change in right end
 * @param ranges : Array of ranges, store isolation results
 * @param num_roots : Store the number of roots
 */
template <int n>
void BudanSquareFreeSolve(const Poly<n> &poly, int duplicate_times, double left,
                          int left_change, double right, int right_change,
                          Range *ranges, int *num_roots) {
  // no root in this range
  if (left_change == right_change)
    return;
  // exact one root in this range
  else if ((right - left) < kMINRANGE && (left_change - right_change) == 1) {
    Range ans{left, right};
    for (int i = 0; i < duplicate_times; i++) {
      ranges[*num_roots] = ans;
      (*num_roots)++;
    }
    return;
  }

  double mid((left + right) / 2.0);
  int mid_change(AddToX(poly, mid).SignChange());

  // left side
  BudanSquareFreeSolve<n>(poly, duplicate_times, left, left_change, mid,
                          mid_change, ranges, num_roots);
  // right side
  BudanSquareFreeSolve<n>(poly, duplicate_times, mid, mid_change, right,
                          right_change, ranges, num_roots);
}

/**
 * @brief Use Budan' Theorem to get range of every root of polynomial that
 * represented by coef and coef_num. Duplicated roots will represented by
 * same range.
 *
 * @param coef :Coefficients of polynomial
 * @param coef_num :Number of coefficients
 * @param ranges :Array of ranges store the isolation results
 * @return :Number of roots
 */
int BudanRootIsolate(const double *coef, int coef_num, Range *ranges) {
  Poly<kMAXDEGREE> original_poly(coef, coef_num);
  Poly<kMAXDEGREE> square_free_polys[kMAXDEGREE];

  int num_roots(0);
  int num_square_free(
      SquareFreeDecompose<kMAXDEGREE>(original_poly, square_free_polys));

  for (int i = 0; i < num_square_free; i++) {

    // Handle zero roots
    HandleZeroRoots(i + 1, &square_free_polys[i], ranges, &num_roots);

    if (square_free_polys[i].get_degree() == 0) // constant
      continue;
    else if (square_free_polys[i].get_degree() == 1) // linear
      HandleLinear<kMAXDEGREE>(square_free_polys[i], i + 1, ranges, &num_roots);
    else if (square_free_polys[i].get_degree() == 2) // quadratic
      HandleQuadratic<kMAXDEGREE>(square_free_polys[i], i + 1, ranges,
                                  &num_roots);
    else { // more than quadratic
      double right(UpperBound(square_free_polys[i])), left = -right;
      int left_change(AddToX(square_free_polys[i], left).SignChange());
      int right_change(AddToX(square_free_polys[i], right).SignChange());

      BudanSquareFreeSolve<kMAXDEGREE>(square_free_polys[i], i + 1, left,
                                       left_change, right, right_change, ranges,
                                       &num_roots);
    }
  }

  return num_roots;
}

#endif // POLY_BUDAN_H
