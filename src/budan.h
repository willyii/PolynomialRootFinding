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
void BudanSquareFreeSolve(const Poly<n> &poly, int duplicate_times,
                          interval left, int left_change, interval right,
                          int right_change, Range *ranges, int *num_roots) {

  // std::cout << "tt search from " << boost::numeric::median(left) << " to "
  //<< boost::numeric::median(right) << std::endl;

  // no root in this range
  if (left_change == right_change)
    return;
  // search range smaller than threshold
  else if ((left_change - right_change) == 1) { // exact one root
    AddToRange(duplicate_times, left, right, ranges, num_roots);
    return;
  } else if (boost::numeric::median(right) - boost::numeric::median(left) <=
             boost::numeric::width(right - left) + 1e-6) {
    // AddToRange(duplicate_times * (left_change - right_change), left, right,
    // ranges, num_roots);
    return;
  }

  interval mid((left + right) / 2.0);
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

  int num_roots(0);
  Poly<kMAXDEGREE> original_poly(coef, coef_num);
  // std::cout << "DEBUG: original poly " << original_poly << std::endl;

  if (original_poly.get_degree() == 1) {
    Linear(original_poly, 1, ranges, &num_roots);
    return 1;
  } else if (original_poly.get_degree() == 2) {
    Quadratic(original_poly, 1, ranges, &num_roots);
    return 2;
  } else {

    Poly<kMAXDEGREE> square_free_polys[kMAXDEGREE];

    int num_square_free(
        SquareFreeDecompose<kMAXDEGREE>(original_poly, square_free_polys));

    for (int i = 0; i < num_square_free; i++) {

      // std::cout << "DEBUG: square free poly " << square_free_polys[i]
      //<< std::endl;

      // Handle zero roots
      if (ZeroRoots(&square_free_polys[i]))
        AddToRange(i + 1, 0.0, 0.0, ranges, &num_roots);

      if (square_free_polys[i].get_degree() == 0) // constant
        continue;
      else if (square_free_polys[i].get_degree() == 1) // linear
        Linear<kMAXDEGREE>(square_free_polys[i], i + 1, ranges, &num_roots);
      else if (square_free_polys[i].get_degree() == 2) // quadratic
        Quadratic<kMAXDEGREE>(square_free_polys[i], i + 1, ranges, &num_roots);
      else {
        interval right(UpperBound(square_free_polys[i])), left = -right;
        int left_change(AddToX(square_free_polys[i], left).SignChange());
        int right_change(AddToX(square_free_polys[i], right).SignChange());

        BudanSquareFreeSolve(square_free_polys[i], i + 1, left, left_change,
                             right, right_change, ranges, &num_roots);
      }
    }

    return num_roots;
  }
}

#endif // POLY_BUDAN_H
