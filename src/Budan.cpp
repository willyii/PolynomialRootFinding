#include "Budan.h"

#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "Param.h"
#include "Util.h"

using namespace std;
using util::coeffAfter;
using util::deg;
using util::isZeroVec;
using util::lc;
using util::polyTimes;

/**
 * Caculate the number of sign change when this equation added h
 *
 * @param poly Polynomila function that need to be caculated
 * @param h after add this number, formula will be changed
 * @return the number of sign change
 */
int Budan::signChangeNums(Poly &poly, double h) {
  int ret = 0;
  vector<double> after = coeffAfter(poly.getCoef(), h);

  bool sign = after[0] > 0;
  for (auto num : after) {
    if (num == 0) continue;
    if (num > 0 == sign) continue;
    sign = !sign;
    ret++;
  }
  return ret;
}

/**
 * NewtonRaphson method to find root.
 *
 * This is already garentee that only one root in polynomoal.
 * TODO: Add intervals to this function and make sure it will find the root.
 *
 * @param poly polynomial function
 * @param x initial guess of the root.
 * @return the root of this polynomial.
 */
double Budan::NewtonRaphson(Poly &poly, double x) {
  int idx = 0;
  while ((abs(poly.getValue(x) - 0) > Param::EPSILON) &&
         (idx < Param::MAX_ITER)) {
    x -= poly.getValue(x) / poly.getGradient(x);
    idx += 1;
  }
  if (abs(poly.getValue(x) - 0) <= Param::EPSILON) return x;
  return numeric_limits<double>::max();
}

/**
 * Finding GCD of two Poly's coefficient by Euclid's algorithm
 *
 * @param c1, c2, coeff of two polys
 * @return gcd of two polys
 */
Poly Budan::gcd(vector<double> &a, vector<double> &b) {
  // Reach the end
  if (isZeroVec(b)) {
    return Poly(a);
  }

  int N = a.size();
  vector<double> q(N, 0.0), r = a;
  int d = deg(b);
  double c = lc(b);

  while (deg(r) >= d) {
    vector<double> s(N, 0.0);
    s[N - deg(r) + d - 1] = lc(r) / c;
    vector<double> sb = polyTimes(s, b);
    for (int i = 0; i < N; i++) {
      q[i] = q[i] + s[i];
      r[i] = r[i] - sb[i];
    }
  }
  return gcd(b, r);
}

/**
 * Apply Yun's Algorithm to decompose a poly in to several square-free poly
 *
 * TODO:
 *
 * @param poly Polynomial need to be decompose
 * @return The coeff of decomposed coef
 */
vector<vector<double>> Budan::squareFreeDecompoe(Poly &poly) {
  vector<Poly> ans;
  Poly a = Budan::gcd(poly.getCoef(), poly.getGradCoef());
  return vector<vector<double>>{{0.0}};
}
