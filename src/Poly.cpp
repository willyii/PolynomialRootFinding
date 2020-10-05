#include "poly.h"

#include <vector>
#include <string.h>

using std::vector;


/**
 * Calculate the coefficient of gradient
 */
void Poly::gradientCoef() {
  int N = _coef.size();
  vector<double> tmp(1,0.0);
  for (int i = 0; i < N - 1; i++) 
    tmp.push_back(_coef[i] * (N - i - 1));
  _gradCoef = Coef(tmp);
  return;
}

// OverLoad Add operator
Poly Poly::operator+(Poly& b) {
  Coef tmp = this->_coef + b.getCoef();
  return Poly(tmp);
}

// OverLoad Sub operator
Poly Poly::operator-(Poly& b) {
  Coef tmp = this->_coef - b.getCoef();
  return Poly(tmp);
}

// OverLoad Times operator
Poly Poly::operator*(Poly& b) {
  Coef tmp = this->_coef * b.getCoef();
  return Poly(tmp);
}

// OverLoad divide operator
Poly Poly::operator/(Poly& b) {
  Coef tmp = this->_coef / b.getCoef();
  return Poly(tmp);
}

// Value at point x
double Poly::valueAt(double x) {
  double ret = 0.0;
  for(int i=0;i<_coef.size();i++)
    ret = ret * x + _coef[i];
  
  return ret;
}

// Gradient at point x
double Poly::gradientAt(double x) {
  double ret = 0.0;
  for(int i=0;i<_gradCoef.size();i++)
    ret = ret * x + _gradCoef[i];
  
  return ret;
}
