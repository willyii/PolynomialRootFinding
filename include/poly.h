#ifndef POLY_H
#define POLY_H

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#include <cassert>
#include <iostream>

template <int n>  // maximum degree, could be less
class Poly {
 private:
  double coef[n + 1];
  int N;  // number of coefficent

 public:
  Poly();
  Poly(double* input_coef, int num_input);

  double valueAt(double x);     // value at point x
  double gradientAt(double x);  // gradient at point x
  int size() const { return N; }

  const double& operator[](int i) const { return coef[i]; }
  template <int n2>
  Poly<n + n2> operator+(Poly<n2>& p);
  // Poly operator-(Poly& p);
  // Poly operator*(Poly& p);
  // Poly operator/(Poly& p);
  // Poly& operator+=(Poly& p);
  // Poly& operator-=(Poly& p);
  // Poly& operator*=(Poly& p);
  // Poly& operator/=(Poly& p);
};

template <int n>
Poly<n>::Poly() {
  for (int i = 0; i < n; i++) coef[i] = 0;
  N = n;
}

template <int n>
Poly<n>::Poly(double* input_coef, int num_input) {
  assert(num_input <= n + 1);
  for (int i = 0; i < num_input; i++) coef[i] = input_coef[i];
  N = num_input;
}

template <int n>
double Poly<n>::valueAt(double x) {
  double ans = 0;
  for (int i = N - 1; i >= 0; i--) ans = ans * x + coef[i];
  return ans;
}

template <int n>
double Poly<n>::gradientAt(double x) {
  double ans = 0;
  for (int i = N - 1; i >= 0; i--) ans = ans * x + coef[i] * i;
  return ans;
}

template <int n, int n2>
Poly<n + n2> Poly<n>::Poly operator+(Poly<n2>& p) {}

template <int n>
static std::ostream& operator<<(std::ostream& out, const Poly<n>& u) {
  for (int i = 0; i < u.size(); i++) {
    out << " ";
    if (u[i] >= 0) out << '+';
    out << u[i] << "x^" << i;
  }
  return out;
}

#endif
