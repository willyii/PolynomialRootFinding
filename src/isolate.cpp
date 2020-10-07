#include "isolate.h"

// GCD of two Poly
Poly Isolate::gcd(Poly& p1, Poly& p2) {
  if (p2.getCoef().isZero()) {
    p1.monic();
    return p1;
  }

  int N = p1.size(), d = p2.getCoef().deg();
  double c = p2.getCoef().lc();
  vector<double> tmp(N, 0.0);
  Poly q = Poly(tmp), r = p1;

  while (r.getCoef().deg() >= d && !r.getCoef().isZero()) {
    Poly s = Poly(tmp);
    s[N - r.getCoef().deg() + d - 1] = r.getCoef().lc() / c;
    Poly sb = s * p2;
    q = q + s;
    r = r - sb;
  }
  return gcd(p2, r);
}

// Square-free decompo
vector<Poly> Isolate::squareFreeDecompo(Poly& p) {
  vector<Poly> ans;
  Poly a, b, c, d, fd, bd;
  fd = p.getGradPoly();
  a = gcd(p, fd);
  b = p / a;
  c = fd / a;
  bd = b.getGradPoly();
  d = c - bd;

  while (!isOne(b)) {
    a = gcd(b, d);
    b = b / a;
    c = d / a;
    bd = b.getGradPoly();
    d = c - bd;
    ans.emplace_back(a);
  }
  return ans;
}

// Check if one
bool Isolate::isOne(Poly& p) {
  int N = p.size();
  for (int i = 0; i < N - 1; i++) {
    if (abs(p[i] - 0) > EPSILON) return false;
  }
  return abs(p[N - 1] - 1.0) <= EPSILON;
}
