#include "Poly.h"

#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include "Util.h"
#include "math.h"
#include "string.h"

using namespace std;
using namespace util;

/**
 * Constructor
 *
 * @param coeff: coefficient of the polynomial
 */
Poly::Poly(vector<double> coeff) {
  _coefficient = coeff;
  GradientCoef();
  _N = coeff.size();
}

/**
 * Calculate the coefficient of gradient
 */
void Poly::GradientCoef() {
  int length = _coefficient.size();
  _gradientCoef = {0.0};
  for (int i = 0; i < length - 1; i++) {
    _gradientCoef.push_back(_coefficient[i] * (length - i - 1));
  }
}

/**
 * Calculate the value of corresponding x
 *
 * @param x: the x value of this equation
 * @return the corresponding value y
 */
double Poly::getValue(double x) {
  double ret = 0.0;
  for (auto tmp : _coefficient) {
    ret = ret * x + tmp;
  }
  return ret;
}

/**
 * Calculate the gradient of poly corresponding x
 *
 * @param x: the x value of this equation
 * @return the corresponding gradient y
 */
double Poly::getGradient(double x) {
  double ret = 0.0;
  for (auto tmp : _gradientCoef) {
    ret = ret * x + tmp;
  }
  return ret;
}

// OverLoad divide operator
Poly Poly::operator/(Poly& b) {
  return Poly(polyDiv(this->getCoef(), b.getCoef()));
}

// OverLoad Add operator
Poly Poly::operator+(Poly& b) {
  return Poly(polyAdd(this->getCoef(), b.getCoef()));
}

// OverLoad Sub operator
Poly Poly::operator-(Poly& b) {
  return Poly(polySub(this->getCoef(), b.getCoef()));
}

// OverLoad Times operator
Poly Poly::operator*(Poly& b) {
  return Poly(polyTimes(this->getCoef(), b.getCoef()));
}

// Format print
void Poly::__str__() {
  int d = deg(this->getCoef());
  string s = "";
  for (int i = this->getN() - d - 1; i < this->getN(); i++) {
    s += " " + to_string(this->getCoef()[i]) + "x^" + to_string(d) + " ";
    d -= 1;
  }
  cout << s << endl;
}