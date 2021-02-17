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

#include "poly.h"
#include "range.h"
#include "util.h"

/**
 * struct Mobius - represent Mobius transformation, map x to (ax + b)/(cx+d)
 */
struct Mobius {
  interval a, b, c, d;
};

/**
 * Help function used in VincentSquareFreeSolve, It will change change
 * p(x) to (x+1)^m p(1/(x+1)), where m is the degree of polynomial
 *
 * @tparam n :Maximum degree of poly
 * @param poly :Polynomial
 * @return :Polynomial after transformation
 */
template <int n> Poly<n> Vincent_helper(const Poly<n> &poly) {
  Poly<n> ret;
  for (int i = poly.get_degree(); i >= 0; i--) {
    ret[poly.get_degree() - i] = poly[i];
  }
  ret.set_degree();
  return AddToX(ret, 1);
}

/**
 * Isolate the roots of square free polynomial using continued fractions
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial
 * @param duplicate_times :Repeat time of this poly repeat in original poly
 * @param mobuis :Mobius transormation to x in poly
 * @param ranges : Store isolation results, might be modified
 * @param num_roots : Store the number of roots, might be modified
 */
template <int n>
void VincentSquareFreeSolve(const Poly<n> &poly, int duplicate_times,
                            const Mobius &mobius, Range *ranges,
                            int *num_roots) {
  int sign_change(poly.SignChange());

  if (sign_change == 0) // no positive root
    return;
  else if (sign_change == 1) { // Only one positive ropot
    AddToRange(duplicate_times, mobius.b / mobius.d, mobius.a / mobius.c,
               ranges, num_roots);
    return;
  }

  auto poly2(AddToX(poly, 1));
  int sign_change_dif(sign_change - poly2.SignChange());

  // Rational root
  if (ZeroRoots(&poly2)) {
    interval root((mobius.a + mobius.b) / (mobius.c + mobius.d));
    AddToRange(duplicate_times, root, root, ranges, num_roots);
    sign_change_dif--;
  }

  // roots in (b, inf)
  Mobius next_mobius{mobius.a, mobius.b + mobius.a, mobius.c,
                     mobius.d + mobius.c};
  VincentSquareFreeSolve<n>(poly2, duplicate_times, next_mobius, ranges,
                            num_roots);

  if (sign_change_dif == 0) // no root in (0,b)
    return;
  else if (sign_change_dif == 1) { // one root in (0,b)
    interval left(mobius.b / mobius.d);
    interval right((mobius.a + mobius.b) / (mobius.c + mobius.d));
    AddToRange(duplicate_times, left, right, ranges, num_roots);
  } else { // more than one roots in  (0,b)

    // change M(x) to M(1/(x+1))
    next_mobius.a = mobius.b;
    next_mobius.b = mobius.a + mobius.b;
    next_mobius.c = mobius.d;
    next_mobius.d = mobius.c + mobius.d;

    VincentSquareFreeSolve<n>(Vincent_helper(poly), duplicate_times,
                              next_mobius, ranges, num_roots);
    return;
  }
}

/**
 * Transform p(x) to p(-x)
 *
 * @tparam n :Maximum degree of polynomial
 * @param poly :Polynomial after transformation
 */
template <int n> Poly<n> NegatePoly(const Poly<n> &poly) {
  Poly<n> ret(poly);
  for (int i = 1; i <= poly.get_degree(); i += 2)
    ret[i] = -1.0 * ret[i];
  return ret;
}

/**
 * Root isolation using continued fractions
 *
 * @param coef :Pointer to coefficients of polynomial
 * @param coef_num :Number of coefficients
 * @param ranges :Store isolation results, might be modified
 * @return :Number of roots
 */
int VincentRootIsolate(const double *coef, int coef_num, Range *ranges) {
  Poly<kMAXDEGREE> original_poly(coef, coef_num);
  // std::cout << "DEBUG: Original polynomial " << original_poly << std::endl;
  Poly<kMAXDEGREE> square_free_polys[kMAXDEGREE];

  // Square free decompose polynomial
  int num_roots(0);
  int num_square_free(
      SquareFreeDecompose<kMAXDEGREE>(original_poly, square_free_polys));

  for (int i = 0; i < num_square_free; i++) {

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

      /* TODO : divied by 0 */

      // Positive roots
      Mobius inital_mobius{1, 0, 0, 1};
      VincentSquareFreeSolve(square_free_polys[i], i + 1, inital_mobius, ranges,
                             &num_roots);
      // Negative roots
      inital_mobius.a = -1;
      VincentSquareFreeSolve(NegatePoly(square_free_polys[i]), i + 1,
                             inital_mobius, ranges, &num_roots);
    }
  }

  return num_roots;
}

#endif // POLY_VINCENT_H
