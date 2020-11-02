#ifndef UTIL_H
#define UTIL_H

#include "poly.h"
#include <unistd.h>
#include <algorithm>
#include "boundry.h"
#include <deque>
#include <tuple>

using std::cout;
using std::deque;
using std::endl;
using std::tuple;
using std::get;

// Compute the GCD between p1 and p2
Poly gcd(Poly p1, Poly p2) {
  if (DEBUG_GCD) {
    std::cout << "DEBUG GCD: \n"
      << "p1: " << p1 << "\t | \t p2: " << p2 << std::endl;
  }
  if (p2.isZero()) {
    return p1;
  }

  p1.monic();
  p2.monic();
  int N = p1.size(), d = p2.deg();
  double c = p2.lc();
  vector<double> tmp(N, 0.0);
  Poly q = Poly(tmp), r = p1;

  while (r.deg() >= d) {
    std::fill(tmp.begin(), tmp.end(), 0.0);
    tmp[N - r.deg() + d - 1] = r.lc() / c;
    Poly s = Poly(tmp);
    Poly sb = s * p2;
    q += s;
    r = r - sb;
  }
  if (DEBUG_GCD) usleep(500000);
  return gcd(p2, r);
}

// Decompose p into square free polynomials
vector<Poly> squareFreeDecompo(Poly& p) {
  if (DEBUG_SQD) std::cout << "DEGUG SQF: p in first: " << p << std::endl;
  vector<Poly> ans;
  Poly a, b, c, d, fd, bd;
  fd = p.getGradPoly();
  if (DEBUG_SQD) std::cout << "DEGUG SQF: fd in first: " << fd << std::endl;
  a = gcd(p, fd);
  if (DEBUG_SQD) std::cout << "DEGUG SQF: a in first: " << a << std::endl;
  b = p / a;
  if (DEBUG_SQD) std::cout << "DEGUG SQF: b in first: " << b << std::endl;
  c = fd / a;
  if (DEBUG_SQD) std::cout << "DEGUG SQF: c in first: " << c << std::endl;
  bd = b.getGradPoly();
  if (DEBUG_SQD) std::cout << "DEGUG SQF: bd in first: " << bd << std::endl;
  d = c - bd;
  if (DEBUG_SQD) std::cout << "DEGUG SQF: d in first: " << d << std::endl;

  while (!b.isOne()) {
    a = gcd(b, d);
    if (DEBUG_SQD) std::cout << "DEGUG SQF: a in others: " << a << std::endl;
    b = b / a;
    if (DEBUG_SQD) std::cout << "DEGUG SQF: b in others: " << b << std::endl;
    c = d / a;
    if (DEBUG_SQD) std::cout << "DEGUG SQF: c in others: " << c << std::endl;
    bd = b.getGradPoly();
    d = c - bd;
    if (DEBUG_SQD) std::cout << "DEGUG SQF: d in others: " << d << std::endl;
    ans.emplace_back(a);
    if (DEBUG_SQD) usleep(300000);
  }
  return ans;
}

// Add convert x to (x + h) return coresponding Poly
Poly addToP(Poly& p, double h) {
  int N = p.size();
  vector<double> tmp(N, 0.0);
  unordered_map<int, Poly> memo;  // map the idx to corresponding coef
  // zero case
  tmp[N - 1] = 1;
  memo[0] = Poly(tmp);
  fill(tmp.begin(), tmp.end(), 0);
  // one case
  if (N >= 1) {
    tmp[N - 1] = h;
    tmp[N - 2] = 1;
    memo[1] = Poly(tmp);
    fill(tmp.begin(), tmp.end(), 0);
  }
  // more case
  for (int idx = 2; idx < N; idx++) {
    int mid = idx / 2;
    memo[idx] = memo[mid] * memo[idx - mid];
  }
  // combine
  Poly ans = Poly(tmp);
  for (int i = 0; i < N; i++) {
    memo[i] = memo[i] * p[N - i - 1];
    ans += memo[i];
  }
  if (DEBUG_BUDAN)
    cout << "DEBUG BUDAN: after add " << h << " : " << Poly(tmp) << endl;
  return ans;
}

// Time h to x, convert it to (hx), return corresponding Polynomial
Poly timeToP(Poly& p, double h) {
  vector<double> coef = p.getCoef();
  double apply = 1.0;
  for(int i= p.size()-1; i>=0 ;i --){
    coef[i] *= apply;
    apply*= h;
  }
  return Poly(coef);
}

// Compute the sign change of coefficient
int signChangeNum(Poly& tmp, double h) {
  Poly p = addToP(tmp, h);
  int ret = 0;
  bool sign = p[0] > 0;
  for (int i = 0; i < p.size(); i++) {
    if (abs(p[i] - 0) <= EPSILON) continue;
    if (p[i] > 0 == sign) continue;
    sign = !sign;
    ret++;
  }
  return ret;
}

// Get the boundry of roots, applied Cauchy's bound
static double bound(Poly& p) {
  int N = p.size();
  double tmp = __DBL_MIN__, lc = p.lc();
  for (int i = 1; i < N; i++) {
    tmp = fmax(tmp, fabs(p[i] / lc));
  }
  return 1 + tmp;
}

// Finding root in boundry b. root is isolated in that boundry
// Newtown method temporary
static double rootInBound(Poly& p, double left, double right) {
  if (DEBUG_BUDAN)
    cout << "DEBUG BUDAN: ===== Searching root in range =====" << left
      << " to " << right << "\n";
  double x0 = (left + right) / 2;
  int idx = 0;
  while (p.valueAt(x0) != 0 && idx < MAXITER && x0 >= left && x0 <= right) {
    double step = (p.valueAt(x0) / p.gradientAt(x0));
    x0 = x0 - step;
    if (step <= EPSILON) break;
    idx += 1;
  }

  if (idx == MAXITER || left - x0 > EPSILON || x0 - right > EPSILON) {
    return NOTFOUND;
  }
  return x0;
}
#endif

