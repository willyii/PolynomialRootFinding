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
      if(p1.deg() <= 0) p1/= p1;
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
      if(DEBUG_GCD) std::cout<<"DEBUG GCD: s = "<<r.lc() / c<<std::endl;
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
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: fd in first: "<<fd<<std::endl;
    a = gcd(p, fd);
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: a in first: "<<a<<std::endl;
    b = p / a;
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: b in first: "<<b<<std::endl;
    c = fd / a;
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: c in first: "<<c<<std::endl;
    bd = b.getGradPoly();
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: bd in first: "<<bd<<std::endl;
    d = c - bd;
    if(DEBUG_SQD) std::cout<<"DEGUG SQF: d in first: "<<d<<std::endl;

    while (!isOne(b)) {
      a = gcd(b, d);
      if(DEBUG_SQD) std::cout<<"DEGUG SQF: a in others: "<<a<<std::endl;
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