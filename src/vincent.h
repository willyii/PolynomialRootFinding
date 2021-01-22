// ----------------------------------------------------------------------------
//
// FILENAME: vincent.h
//
// DESCRIPTION:
//    This file define and implement functions that will be used in root
//    isolation with continued fractions.
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_VINCENT_H
#define POLY_VINCENT_H

#include "range.h"
#include "util.h"

/**
 * @brief :Isolate the roots of square free polynomial using continued fractions
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param duplicate_times :How many times thie poly repeat in original poly
 * @param ranges : Store isolation results, might be modified
 * @param num_roots : Store the number of roots, might be modified
 */
template <int n>
void VincentSquareFreeSolve(const Poly<n> &poly, int duplicate_times,
                            Range *ranges, int *num_roots) {
  int sign_change(poly.SignChange());
  /* TODO : */
}

/**
 * @brief :Isolate roots using continued fractions
 *
 * @param coef :Coefficent of polynomial
 * @param coef_num :Number of coefficients
 * @param ranges :Array of range to store result, might be modified
 * @return :Number of roots
 */
int VincentRootIsolate(const double *coef, int coef_num, Range *ranges) {
  Poly<kMAXDEGREE> original_poly(coef, coef_num);
  Poly<kMAXDEGREE> square_free_polys[kMAXDEGREE];

  // Square free decompose polynomial
  int num_roots(0);
  int num_square_free(
      SquareFreeDecompose<kMAXDEGREE>(original_poly, square_free_polys));

  for (int i = 0; i < num_square_free; i++) {
    std::cout << "Debug: " << square_free_polys[i] << std::endl;

    // Handle zero roots
    HandleZeroRoots(i + 1, &square_free_polys[i], ranges, &num_roots);

    if (square_free_polys[i].get_degree() == 0) // constant
      continue;
    else if (square_free_polys[i].get_degree() == 1) // linear
      HandleLinear<kMAXDEGREE>(square_free_polys[i], i + 1, ranges, &num_roots);
    else if (square_free_polys[i].get_degree() == 2) // quadratic
      HandleQuadratic<kMAXDEGREE>(square_free_polys[i], i + 1, ranges,
                                  &num_roots);
    else {
      /* TODO :  */
    }
  }

  return num_roots;
}

#endif // POLY_VINCENT_H
