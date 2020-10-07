#include "poly.h"

#include <string.h>

#include <vector>

using namespace std;

// Constructor
Poly::Poly() {}
Poly::Poly(Coef& d) {
  _coef = d;
  gradientCoef();
}
Poly::Poly(vector<double>& d) {
  _coef = Coef(d);
  gradientCoef();
}

// Accessor
Coef& Poly::getCoef() { return _coef; }
Coef& Poly::getGradCoef() { return _gradCoef; }
int Poly::size() const { return _coef.size(); }
Poly Poly::getGradPoly() { return Poly(_gradCoef); }

// Overload [] operator
const double& Poly::operator[](int i) const { return _coef[i]; }
double& Poly::operator[](int i) { return _coef[i]; }

// Calculate the coefficient of gradient
void Poly::gradientCoef() {
  int N = _coef.size();
  vector<double> tmp(1, 0.0);
  for (int i = 0; i < N - 1; i++) tmp.push_back(_coef[i] * (N - i - 1));
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
  for (int i = 0; i < _coef.size(); i++) ret = ret * x + _coef[i];

  return ret;
}

// Gradient at point x
double Poly::gradientAt(double x) {
  double ret = 0.0;
  for (int i = 0; i < _gradCoef.size(); i++) ret = ret * x + _gradCoef[i];

  return ret;
}

// Overload << operator
std::ostream& operator<<(std::ostream& out, const Poly& u) {
  int N = u.size();
  for (int i = 0; i < N; i++) {
    if (i) out << ' ' << '+' << ' ';
    out << u[i] << 'x' << '^' << N - i - 1;
  }
  return out;
}

// Make polynomial monic
void Poly::monic() {
  double div = 1;
  for (int i = 0; i < _coef.size(); i++) {
    if (fabs(_coef[i] - 0.0) > EPSILON) {
      div = 1 / _coef[i];
      break;
    }
  }
  for (int i = 0; i < _coef.size(); i++) {
    if (fabs(_coef[i] - 0.0) > EPSILON) {
      _coef[i] *= div;
    }
  }
  return;
}