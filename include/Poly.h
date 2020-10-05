#ifndef POLY_H
#define POLY_H

#include "coef.h"
#include <iostream>

using std::vector;

class Poly {
  public:
  Poly() {};
  Poly(vector<double>& d) {_coef = Coef(d); gradientCoef();} 
  Poly(Coef& d){_coef = d; gradientCoef();}
  Coef& getCoef() { return _coef; }
  Coef& getGradCoef() { return _gradCoef; }
  int size() const {return _coef.size();}
  double valueAt(double x);
  double gradientAt(double x);

  Poly operator+(Poly& b);
  Poly operator-(Poly& b);
  Poly operator*(Poly& b);
  Poly operator/(Poly& b);
  const double& operator[](int i) const {return _coef[i];}

  Poly getGradPoly() {return Poly(_gradCoef);}


  private:
  void gradientCoef();
  coef _coef;
  coef _gradCoef;


};

// Overload << operator
static std::ostream& operator<<(std::ostream& out, const Poly& u){
  int N = u.size();
  for (int i = 0; i < N; i++) {
    if (i) out << ' ' << '+' << ' ';
    out << u[i] << 'x' <<'^' << N-i-1;
  }
  return out;
}

#endif
