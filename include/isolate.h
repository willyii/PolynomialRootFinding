#ifndef ISOLATE_H
#define ISOLATE_H

#include "poly.h"

using std::vector;

class Isolate {
 public:
  static Poly gcd(Poly& p1, Poly& p2);
  static vector<Poly> squareFreeDecompo(Poly& p);

 protected:
  static bool isOne(Poly& p);
};

#endif