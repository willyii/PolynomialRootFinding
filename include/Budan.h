#ifndef BUDAN_H
#define BUDAN_H

#include <vector>

#include "Poly.h"

class Budan {
 public:
  int signChangeNums(Poly &poly, double h);
  double NewtonRaphson(Poly &poly, double x);
  vector<double> rootFinding(Poly &poly);

 private:
  vector<double> coeffAfter(Poly &poly, double h);
  vector<double> polyTimes(vector<double> &c1, vector<double> &c2);
};

#endif