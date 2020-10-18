#ifndef ISOLATE_H
#define ISOLATE_H

#include "poly.h"

class Isolate {
 public:
  static Poly gcd(Poly& p1, Poly& p2) {
    if (DEBUG_GCD) {
      std::cout << "DEBUG GCD: \n"
                << "p1: " << p1 << "\t | \t p2: " << p2 << std::endl;
    }
    if (p2.isZero()) {
      return p1;
    }

    int N = p1.size(), d = p2.deg();
    double c = p2.lc();
    vector<double> tmp(N, 0.0);
    Poly q = Poly(tmp), r = p1;

    while (r.deg() >= d && !r.isZero()) {
      std::fill(tmp.begin(), tmp.end(), 0.0);
      tmp[N - r.deg() + d - 1] = r.lc() / c;
      Poly s = Poly(tmp);
      Poly sb = s * p2;
      q += s;
      r = r - sb;
    }
    return gcd(p2, r);
  }
  static vector<Poly> squareFreeDecompo(Poly& p) {
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

 protected:
  static bool isOne(Poly& p) {
    int N = p.size();
    for (int i = 0; i < N - 1; i++) {
      if (fabs(p[i] - 0) > EPSILON) return false;
    }
    return fabs(p[N - 1] - 1.0) <= EPSILON;
  }
};

#endif