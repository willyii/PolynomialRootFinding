#ifndef POLY_H
#define POLY_H

#include <iostream>

#include "coef.h"

using std::vector;

class Poly {
 public:
  Poly();
  Poly(vector<double>& d);
  Poly(Coef& d);
  Coef& getCoef();
  Coef& getGradCoef();
  int size() const;
  double valueAt(double x);
  double gradientAt(double x);

  Poly operator+(Poly& b);
  Poly operator-(Poly& b);
  Poly operator*(Poly& b);
  Poly operator/(Poly& b);
  friend std::ostream& operator<<(std::ostream& out, const Poly& u);
  const double& operator[](int i) const;
  double& operator[](int i);

  Poly getGradPoly();
  void monic();

 private:
  void gradientCoef();
  coef _coef;
  coef _gradCoef;
};

#endif
