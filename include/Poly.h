#ifndef POLY_H
#define POLY_H

#include "coef.h"
#include <iostream>

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

  Poly getGradPoly();


  private:
  void gradientCoef();
  coef _coef;
  coef _gradCoef;


};

#endif
