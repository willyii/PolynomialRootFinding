#ifndef UTIL_H
#define UTIL_H

//#include <unistd.h>

//#include <algorithm>
//#include <deque>
//#include <tuple>

#include "poly.h"

// Compute the GCD between p1 and p2
// /* TODO: */
// Poly gcd(Poly p1, Poly p2);

// Decompose p into square free polynomials
// /* TODO: */
// vector<Poly> squareFreeDecompo(Poly& p);

// Add convert x to (x + h) return coresponding Poly
// /* TODO: */
template <int n>
Poly<n> addToX(Poly<n>& p, double h) {
  Poly<n> ret;

  Poly<0> op(1.0, 0);

  Poly<1> xh(h, 0);
  xh.addCoef(1, 1);

  for (int i = 0; i < p.size(); i++) {
    ret = ret + op * p[i];
    // op = op * xh;
  }

  return ret;
}

// Time h to x, convert it to (hx), return corresponding Polynomial
template <int n>
Poly<n> timeToX(Poly<n>& p, double h) {
  Poly<n> ret(p);
  double op = 1;
  for (int i = 0; i < ret.size(); i++) {
    ret[i] *= op;
    op *= h;
  }
  return ret;
}

// Compute the sign change of coefficient
template <int n>
int signChangeNum(Poly<n>& p) {
  int ans = 0;
  bool sign = p[p.degree()] > 0;
  for (int i = p.degree() - 1; i >= 0; i--)
    if (fabs(p[i]) > EPSILON && ((p[i] > 0) != sign)) {
      sign = !sign;
      ans++;
    }
  return ans;
}

// Get the upper boundry of roots, applied Cauchy's bound
template <int n>
double upperBound(Poly<n>& p) {
  double ans = 0, lc = p.lc();
  for (int i = 0; i < p.degree(); i++) ans = fmax(ans, fabs(p[i] / lc));
  return 1 + ans;
}

// Get the lower boundary of roots, applied Cauchy's bound
template <int n>
double lowerBound(Poly<n>& p) {
  double a0 = p[0], ans = 0;
  if (fabs(a0) < EPSILON) return 0.0;
  for (int i = 1; i < p.size(); i++) ans = fmax(ans, fabs(p[i] / a0));
  return 1 / (1 + ans);
}

// Finding root in boundry b. root is isolated in that boundry
// Newtown method temporary
/* TODO: */
// double rootInBound(Poly& p, double left, double right);

#endif
