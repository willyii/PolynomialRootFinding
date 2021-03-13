// ----------------------------------------------------------------------------
// FILENAME: timetest.cpp
//
// DESCRIPTION:
//    This file contians the function that used for save the running time of
//    computing polynomial real root isolation
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include "budan.h"
#include "poly.h"
#include "range.h"
#include "vincent.h"
#include <boost/numeric/interval/utility_fwd.hpp>
#include <chrono>
#include <time.h>

static const int kTESTDEGREE = 7;
static const int digit = 4; // number of digit after point
static const int digit_control = std::pow(10, digit); // controler of digit

static const double max_root = 10000;
// static const double max_root = std::pow(2, kTESTDEGREE);

/**
 * Get random double in range min to max
 */
double rand_double(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  f = min + f * (max - min);
  f = std::ceil(f * digit_control) / digit_control;
  return f;
}

int main() {
  srand(time(NULL));

  // Get random polynomial
  double *coeffs = new double[kTESTDEGREE];

  double budan_total = 0, vincent_total = 0;

  for (int i = 0; i < 1000; i++) {

    for (size_t i = 0; i < kTESTDEGREE; i++) {
      coeffs[i] = rand_double(-max_root, max_root);
    }

    Poly<kTESTDEGREE + 1> tt(coeffs, kTESTDEGREE);

    std::cout << tt << std::endl;

    // save roots
    Range *roots = new Range[kTESTDEGREE];

    // Budan
    auto budan_start = std::chrono::high_resolution_clock::now();
    BudanRootIsolate(coeffs, kTESTDEGREE, roots);
    auto budan_end = std::chrono::high_resolution_clock::now();

    // Vincent
    auto vincent_start = std::chrono::high_resolution_clock::now();
    VincentRootIsolate(coeffs, kTESTDEGREE, roots);
    auto vincent_end = std::chrono::high_resolution_clock::now();

    //     Time
    auto budan_duration = std::chrono::duration_cast<std::chrono::microseconds>(
        budan_end - budan_start);
    auto vincent_duration =
        std::chrono::duration_cast<std::chrono::microseconds>(vincent_end -
                                                              vincent_start);

    budan_total += budan_duration.count();
    vincent_total += vincent_duration.count();
  }

  std::cout << "Budan Theorem takes " << budan_total / 1000.0 << " us for "
            << kTESTDEGREE - 1 << " degree" << std::endl;

  std::cout << "Continued Fraction takes " << vincent_total / 1000.0
            << " us for " << kTESTDEGREE - 1 << " degree" << std::endl;

  return 0;
}
