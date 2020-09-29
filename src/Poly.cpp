#include "Poly.h"

#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include "math.h"

using namespace std;

/**
 * Constructor
 *
 * @param coeff: coefficient of the polynomial
 */
Poly::Poly(vector<double> coeff) {
  _coefficient = coeff;
  GradientCoef();
  _N = coeff.size();
}

/**
 * Calculate the coefficient of gradient
 */
void Poly::GradientCoef() {
  int length = _coefficient.size();
  _gradientCoef = {0.0};
  for (int i = 0; i < length - 1; i++) {
    _gradientCoef.push_back(_coefficient[i] * (length - i - 1));
  }
}

/**
 * Calculate the value of corresponding x
 *
 * @param x: the x value of this equation
 * @return the corresponding value y
 */
double Poly::getValue(double x) {
  double ret = 0.0;
  for (auto tmp : _coefficient) {
    ret = ret * x + tmp;
  }
  return ret;
}

/**
 * Calculate the gradient of poly corresponding x
 *
 * @param x: the x value of this equation
 * @return the corresponding gradient y
 */
double Poly::getGradient(double x) {
  double ret = 0.0;
  for (auto tmp : _gradientCoef) {
    ret = ret * x + tmp;
  }
  return ret;
}
