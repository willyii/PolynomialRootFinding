// ----------------------------------------------------------------------------
//
// FILENAME: randompoly.h
//
// DESCRIPTION:
//    This file defines and implement functions used to generate random
//    polynomial with no repeat roots. This polynomials will be writen in to
//    path "kRANDOM_FILE" defined in param.h
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>

#include "param.h"
#include "poly.h"

/**
 * Get random integer in range min to max
 */
int rand_int(int min, int max) { return rand() % (max - min) + min; }

/**
 * Get random double in range min to max
 */
double rand_double(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  return min + f * (max - min);
}

/**
 * Generate a random polynomial.
 *
 */
Poly<kMAXDEGREE> RandomPoly() {
  /* TODO : Generate polynomial with specific root */
  int degree = rand_int(3, kMAXDEGREE - 1);
  double coef[kMAXDEGREE + 1];

  for (int i = 0; i <= degree; i++)
    coef[i] = rand_double(-50, 50);

  Poly<kMAXDEGREE> ans(coef, degree + 1);

  return ans;
}

/**
 * Return corresponding string of polynomial
 */
template <int n> std::string PolyToString(const Poly<n> &poly) {
  std::string ret = "";
  for (int i = 0; i <= n; i++) {
    // if (std::fabs(poly[i]) < kEPSILON)
    //  continue;
    if (poly[i] >= 0)
      ret += '+';
    ret += std::to_string(poly[i]) + "*x^" + std::to_string(i);
  }
  return ret;
}

/**
 * Generate some polynomials and save to file
 *
 * @param num_polys :Number of polynomials to generate
 */
void RandomPolyToFile(int num_polys = 100000) {
  srand(time(NULL));

  std::ofstream write_file;
  write_file.open(kRANDOM_FILE);

  for (int i = 0; i < num_polys; i++) {
    auto single_poly = RandomPoly();
    write_file << PolyToString(single_poly) << "\n";
  }

  write_file.close();
  return;
}
