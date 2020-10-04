#ifndef BUDAN_H
#define BUDAN_H

#include <vector>

#include "Poly.h"


class Budan {
 public:
  int signChangeNums(Poly &poly, double h);
  double NewtonRaphson(Poly &poly, double x);
  Poly gcd(vector<double> &c1, vector<double> &c2);
  vector<vector<double>> squareFreeDecompoe(Poly &poly);

};

#endif
