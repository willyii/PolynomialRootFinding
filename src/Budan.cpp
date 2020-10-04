#include "Budan.h"

#include <unistd.h>

#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "Param.h"
#include "Util.h"
using namespace std;
using namespace util;

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

  while (deg(r) >= d && !isZeroVec(r)) {
    sleep(1);
    vector<double> s(N, 0.0);
    s[N - deg(r) + d - 1] = lc(r) / c;
    vector<double> sb = polyTimes(s, b);
    q = polyAdd(q, s);
    r = polySub(r, sb);
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
vector<Poly> Budan::squareFreeDecompoe(Poly &poly) {
  vector<Poly> ans;
  Poly a, b, c, d;
  a = gcd(poly.getCoef(), poly.getGradCoef());         // a0 = gcd(f, f')
  b = poly / a;                                        // b1 = f/a0
  c = Poly(polyDiv(poly.getGradCoef(), a.getCoef()));  // c = f'/a0
  d = Poly(polySub(c.getCoef(), b.getGradCoef()));     // d = c -b'

  while (!(isOne(b.getCoef()))) {
    a = gcd(b.getCoef(), d.getCoef());
    b = b / a;
    c = Poly(polyDiv(d.getCoef(), a.getCoef()));
    d = Poly(polySub(c.getCoef(), b.getGradCoef()));
    ans.emplace_back(a);
  }
  return ans;
}
