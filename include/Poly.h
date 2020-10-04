#ifndef POLY_H
#define POLY_H

#include <vector>

using std::vector;

class Poly {
 public:
  // Constructor, mutator, accessor
  Poly(){};
  Poly(vector<double> coeff);
  vector<double>& getCoef() { return _coefficient; }
  vector<double>& getGradCoef() { return _gradientCoef; }
  int getN() { return _N; }
  Poly operator/(Poly& b);
  Poly operator+(Poly& b);
  Poly operator-(Poly& b);
  Poly operator*(Poly& b);

  void __str__();

  // Get corresponding value and gradient
  double getValue(double x);
  double getGradient(double x);

 private:
  void GradientCoef();
  
  // variables
  vector<double> _coefficient;
  vector<double> _gradientCoef;
  int _N;  // length of coeff
};

#endif
