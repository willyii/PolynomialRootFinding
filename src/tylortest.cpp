// ----------------------------------------------------------------------------
// FILENAME: tylortest.cpp
//
// DESCRIPTION:
//    This file contians the function that used for save the running time of
//    tylor expansion of polynomial shift with original one
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#include "poly.h"
#include <boost/numeric/interval/utility_fwd.hpp>
#include <chrono>
#include <time.h>

static const int kTYLORDEGREE = 6;
static const int digit = 2; // number of digit after point
static const int digit_control = std::pow(10, digit); // controler of digit
static const double max_root = 1000;
// static const double max_root = std::pow(2, kTYLORDEGREE);

/**
 * Get random double in range min to max
 */
double rand_double(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  f = min + f * (max - min);
  f = std::ceil(f * digit_control) / digit_control;
  return f;
}

/**
 * Replace "x" in polynomial with "x+h"
 * Applied Taylor Expansion to this.
 * p(x+h) = p(h) + p'(h)x + 1/2*p''(h)x^2 ... 1/(n!) * p^n(h)*x^n
 * ref: https://math.stackexchange.com/questions/694565/polynomial-shift
 *
 * @tparam n :Maximum degree of poly
 * @param poly :Polynomial
 * @param h :Number add to x
 * @return :Polynomial that replace x in poly to x+h
 */
template <int n> Poly<n> Tylor(const Poly<n> &poly, interval h) {
  if (h.lower() <= 0.0 && h.upper() >= 0.0)
    return poly;

  Poly<n> ret, tmp(poly);
  ret[0] = tmp.ValueAt(h);
  double divisor(1.0);
  for (int index = 1; index <= poly.get_degree(); ++index) {
    divisor *= index;
    tmp.Derivative_();
    ret[index] = (tmp.ValueAt(h) / divisor);
  }

  ret.set_degree();
  return ret;
}

/**
 * Replace "x" in polynomial with "x+h"
 * With original implementation. Start with x+h, then multiply them
 *
 * @tparam n :Maximum degree of poly
 * @param poly :Polynomial
 * @param h :Number add to x
 * @return :Polynomial that replace x in poly to x+h
 */
template <int n> Poly<n> Original(const Poly<n> &poly, interval h) {
  double tmp[2] = {boost::numeric::median(h), 1};
  Poly<n> tmp_poly(tmp, 2), ret, multiplier(tmp, 2);
  ret[0] = poly[0];

  for (int i = 1; i <= poly.get_degree(); i++) {
    ret += (poly[i] * multiplier);
    auto intermid = multiplier;
    for (int j = multiplier.get_degree(); j >= 0; j--) {
      multiplier[j + 1] = multiplier[j];
    }
    multiplier[0] = 0;
    multiplier += h * intermid;
    multiplier.set_degree();
  }

  ret.set_degree();

  return ret;
}

int main() {
  srand(time(NULL));

  double *coeffs = new double[kTYLORDEGREE];
  for (size_t i = 0; i < kTYLORDEGREE; i++) {
    coeffs[i] = rand_double(-kTYLORDEGREE, kTYLORDEGREE);
  }

  interval h(rand_double(-max_root, max_root));

  Poly<kTYLORDEGREE> test_poly(coeffs, kTYLORDEGREE);
  std::cout << "orig " << test_poly << std::endl;

  // Tylor
  auto tylor_start = std::chrono::high_resolution_clock::now();
  auto r1 = Tylor(test_poly, h);
  auto tylor_end = std::chrono::high_resolution_clock::now();

  //    Original
  auto ori_start = std::chrono::high_resolution_clock::now();
  auto r2 = Original(test_poly, h);
  auto ori_end = std::chrono::high_resolution_clock::now();

  //     Time
  auto tylor_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      tylor_end - tylor_start);
  auto ori_duration = std::chrono::duration_cast<std::chrono::microseconds>(
      ori_end - ori_start);

  std::cout << "Tylor expansion method takes " << tylor_duration.count()
            << " ns for " << kTYLORDEGREE - 1 << " degree" << std::endl;

  std::cout << "Original method takes " << ori_duration.count() << " ns for "
            << kTYLORDEGREE - 1 << " degree" << std::endl;

  std::cout << r1 << std::endl;
  std::cout << r2 << std::endl;

  interval a = 1.0, b = 3.0;
  interval tmp = a - b;
  std::cout << boost::numeric::median(tmp) << "[" << boost::numeric::width(tmp)
            << "]" << std::endl;

  return 0;
}
