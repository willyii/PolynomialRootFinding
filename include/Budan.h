#ifndef BUDAN_H
#define BUDAN_H

#include <vector>

#include "Poly.h"

struct GCD{
  vector<double> q;
  vector<double> r;
};

class Budan {
 public:
  int signChangeNums(Poly &poly, double h);
  double NewtonRaphson(Poly &poly, double x);
  vector<double> rootFinding(Poly &poly);
  GCD gcd(vector<double> &c1, vector<double> &c2);

 private:
  double lc(vector<double> &c);
  int deg(vector<double> &c);
  vector<double> coeffAfter(Poly &poly, double h);
  vector<double> polyTimes(vector<double> &c1, vector<double> &c2);
};

#endif
