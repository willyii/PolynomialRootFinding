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
  getGradientCoef();
  _N = coeff.size();
}

/**
 * Calculate the coefficient of gradient
 */
void Poly::getGradientCoef() {
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

/**
 * NewtonRaphson method to find root.
 *
 * This is already garentee that only one root in polynomoal.
 * TODO: Add intervals to this function and make sure it will find the root.
 *
 * @param x initial guess of the root.
 * @return the root of this polynomial.
 */
double Poly::NewtonRaphson(double x) {
  int idx = 0;
  while ((abs(getValue(x) - 0) > Param::EPSILON) && (idx < Param::MAX_ITER)) {
    x -= getValue(x) / getGradient(x);
    idx += 1;
  }
  if (abs(getValue(x) - 0) <= Param::EPSILON) return x;
  return numeric_limits<double>::max();
}

/**
 * Caculate the number of sign change when this equation added h
 *
 * @param h after add this number, formula will be changed
 * @return the number of sign change
 */
int Poly::signChangeNums(double h) {
  int ret = 0;
  vector<double> after = coeffAfter(h);
  bool sign = after[0] > 0;
  for (auto num : after) {
    if (num == 0) continue;
    if (num > 0 && sign) continue;
    sign = !sign;
    ret++;
  }
  return ret;
}

/**
 * Caculate the coefficient after add h to x
 *
 * @param h the change of x
 * @return coef after change
 */
vector<double> Poly::coeffAfter(double h) {
  vector<double> tmp(_N, 0);
  unordered_map<int, vector<double>> memo;  // map the idx to corresponding coef

  // zero case
  memo[0] = vector<double>(_N, 0);
  memo[0][_N - 1] = 1;

  // one case
  if (_N >= 1) {
    memo[1] = vector<double>(_N, 0);
    memo[1][_N - 1] = h;
    memo[1][_N - 2] = 1;
  }

  // more case
  for (int idx = 2; idx < _N; idx++) {
    int mid = idx / 2;
    memo[idx] = polyTimes(memo[mid], memo[idx - mid]);
  }

  // combine
  for (int i = 0; i < _N; i++) {
    for (int j = 0; j < _N; j++) {
      tmp[j] += _coefficient[_N - i - 1] * memo[i][j];
    }
  }

  return tmp;
}

/**
 * Calculate the coef after two poly times
 * We can make sure the maximum will no large than _N
 *
 * @param c1 coef of poly1
 * @param c2 coef of poly2
 * @return the cofe after computation
 */
vector<double> Poly::polyTimes(vector<double> &c1, vector<double> &c2) {
  vector<double> ret(_N, 0);
  int idx1, idx2;
  for (int i = 0; i < _N; i++) {
    idx1 = _N - i - 1;
    for (int j = 0; j < _N; j++) {
      idx2 = _N - j - 1;
      if (idx1 + idx2 >= _N) continue;
      ret[_N - 1 - idx1 - idx2] += c1[i] * c2[j];
    }
  }
  return ret;
}