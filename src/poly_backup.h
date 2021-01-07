#ifndef POLY_H
#define POLY_H

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#include <cassert>
#include <iostream>

#include "param.h"

template <int n>  // maximum degree, could be less
class Poly {
 private:
  double coef[n + 1];
  int N;  // number of coefficents

 public:
  Poly();
  Poly(double* input_coef, int num_input);
  Poly(double num, int idx);

  template <int n1>
  Poly(const Poly<n1>& p) : coef{} {
    static_assert(n >= n1);
    N = p.size();
    for (int i = 0; i < N; i++) coef[i] = p[i];
  };

  double valueAt(double x);     // value at point x
  double gradientAt(double x);  // gradient at point x
  int size() const { return N; }
  int degree() const { return N - 1; }
  double lc() const { return coef[N - 1]; }

  void setN(int N_in) { N = N_in; }
  void setN() {
    while (N > 1 && fabs(coef[N - 1]) < EPSILON) {
      N--;
    }
  }

  /* TODO:  Get rid of this  */
  void addCoef(int idx, double num) {
    static_assert(idx <= n + 1);
    coef[idx] += num;
    N = max(idx + 1, N);
  }
  void rightShift(int unit, double num) {
    for (int i = N - 1; i >= 0; i--) coef[i + unit] = coef[i] * num;
    for (int i = unit - 1; i >= 0; i--) coef[i] = 0.0;
    N += unit;
  }

  const double& operator[](int i) const { return coef[i]; }
  double& operator[](int i) { return coef[i]; }
  /* TODO: if constexpr(n1>=n2) return Poly<n1>(p1)+=p2; else return
   * Poly<n2>(p2)+=p1; */
  template <int n2>
  Poly<n>& operator+=(Poly<n2> p) {}
};

template <int n>
Poly<n>::Poly() : coef{} {
  N = 1;
}

template <int n>
Poly<n>::Poly(double num, int idx) : coef{} {
  coef[idx] = num;
}

template <int n>
Poly<n>::Poly(double* input_coef, int num_input) : coef{} {
  assert(num_input <= n + 1);
  for (int i = 0; i < num_input; i++) coef[i] = input_coef[i];
  N = num_input;
}

template <int n>
double Poly<n>::valueAt(double x) {
  double ans = coef[N - 1];
  for (int i = N - 2; i >= 0; i--) ans = ans * x + coef[i];
  return ans;
}

template <int n>
double Poly<n>::gradientAt(double x) {
  double ans = coef[N - 1] * (N - 1);
  for (int i = N - 2; i >= 1; i--) ans = ans * x + coef[i] * i;
  return ans;
}

/* TODO: implement += return Poly<...>(a)+=b; */
template <int n1, int n2>
Poly<max(n1, n2)> operator+(const Poly<n1>& p1, const Poly<n2>& p2) {
  Poly<max(n1, n2)> ret(p1);
  for (int i = 0; i < p2.size(); i++) ret[i] += p2[i];

  /*  Set N   */
  if (p1.size() != p2.size())
    ret.setN(max(p1.size(), p2.size()));
  else
    ret.setN();
  return ret;
}

template <int n1, int n2>
Poly<max(n1, n2)> operator-(const Poly<n1>& p1, const Poly<n2>& p2) {
  Poly<max(n1, n2)> ret(p1);
  for (int i = 0; i < p2.size(); i++) ret[i] -= p2[i];

  /*  Set N   */
  if (p1.size() != p2.size())
    ret.setN(max(p1.size(), p2.size()));
  else
    ret.setN();
  return ret;
}

template <int n1, int n2>
Poly<n1 + n2> operator*(const Poly<n1>& p1, const Poly<n2>& p2) {
  Poly<n1 + n2> ret;
  for (int i = 0; i < p1.size(); i++)
    for (int j = 0; j < p2.size(); j++) ret[i + j] += p1[i] * p2[j];
  ret.setN(p1.size() + p2.size() - 1);
  return ret;
}

template <int n1>
Poly<n1> operator*(const Poly<n1>& p1, double num) {
  Poly<n1> ret(p1);
  for (int i = 0; i < p1.size(); i++) ret[i] *= num;
  return ret;
}

template <int n>
struct divRet {
  Poly<n> q;
  Poly<n> r;
};

template <int n1, int n2>
divRet<n1> operator/(const Poly<n1>& p1, const Poly<n2>& p2) {
  Poly<n1> q;
  Poly<n1> r(p1);
  int d = p2.degree(), dr = r.degree();
  double c = p2.lc();
  while (dr >= d) {
    double num = r.lc() / c;
    int idx = dr - d;
    q.addCoef(idx, num);

    Poly<n1> tmp(p2);
    tmp.rightShift(idx, num);
    r = r - tmp;
    dr = r.degree();
  }
  return divRet<n1>{q, r};
}

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
