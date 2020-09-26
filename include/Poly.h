#ifndef POLY_H
#define POLY_H

#include <vector>

#include "Param.h"
using std::vector;

class Poly {
 public:
  Poly(vector<double> coeff);
  Poly(){};
  double getValue(double x);
  double getGradient(double x);
  double NewtonRaphson(double x0);
  int signChangeNums(double h);

 private:
  void getGradientCoef();
  vector<double> coeffAfter(double h);
  vector<double> polyTimes(vector<double> &c1, vector<double> &c2);

  vector<double> _coefficient;
  vector<double> _gradientCoef;
  int _N;
};

#endif