#ifndef POLY_H
#define POLY_H

#include <iostream>
#include <unordered_map>
#include <vector>

#include "math.h"
#include "param.h"

using std::unordered_map;
using std::vector;

class Poly {
 private:
  vector<double> _coeff;
  vector<double> _gradient;
  vector<double> _b;

  void _gradCalc() {
    int N = _coeff.size();
    _gradient = {};
    for (int i = 0; i < N - 1; i++) {
      _gradient.emplace_back((N - i - 1) * _coeff[i]);
    }
  }
  void _cleanUp() {
    int N = 0;
    for (int i = 0; i < _coeff.size(); i++) {
      // if (fabs(_coeff[i] - 0) > EPSILON) break;
      if (_coeff[i] != 0) break;
      N++;
    }
    if (N == _coeff.size()) {
      _coeff = {0.0};
      _gradient = {0.0};
      return;
    }
    _coeff.erase(_coeff.begin(), _coeff.begin() + N);
    _gradient.erase(_gradient.begin(), _gradient.begin() + N);
    return;
  }

 public:
  Poly(){};
  Poly(vector<double>& d) : _coeff(d) {
    _gradCalc();
    _cleanUp();
  }

  void addB(double b) { _b.emplace_back(b); }

  // Basic properties
  vector<double>& getCoef() { return _coeff; }
  vector<double>& getGradCoef() { return _gradient; }
  int size() const { return _coeff.size(); }
  bool isZero() const {
    int N = size();
    for (int i = 0; i < N; i++)
      if (fabs(_coeff[i] - 0.0) > EPSILON) return false;
    return true;
  }
  bool isOne() const {
    int N = size();
    for (size_t i = 0; i < N - 1; i++)
      if (fabs(_coeff[i] - 0.0) > EPSILON) return false;
    return fabs(_coeff[N - 1] - 1.0) <= EPSILON;
  }

  int deg() const {
    int N = _coeff.size(), i;
    for (i = 0; i < N; i++) {
      if (fabs(_coeff[i] - 0.0) > EPSILON) return N - 1 - i;
    }
    return -1;
  }
  double lc() const {
    for (int i = 0; i < size(); i++)
      if (fabs(_coeff[i] - 0.0) > EPSILON) return _coeff[i];
    return 0;
  }
  Poly leadDiv(Poly& c2) {
    double div = lc() / c2.lc();
    int degree = deg() - c2.deg();
    vector<double> t(size(), 0.0);
    t[size() - degree - 1] = div;
    return Poly(t);
  }
  Poly getGradPoly() { return Poly(_gradient); }
  void monic() {
    double tmp = lc();
    for (int i = 0; i < _coeff.size(); i++) {
      _coeff[i] /= tmp;
    }
    _gradCalc();
    _cleanUp();
  };
  Poly pow(const int d) {
    unordered_map<int, Poly> memo;
    vector<double> tmpCoef = {1.0};
    memo[0] = Poly(tmpCoef);
    memo[1] = *this;
    for (size_t i = 2; i <= d; i++) {
      int mid = i / 2;
      memo[i] = memo[mid] * memo[i - mid];
    }
    return memo[d];
  }

  // Some calculation
  double valueAt(double x) {
    if (_b.size() > 0) {
      for (auto b : _b) {
        x = b / (1 + x);
      }
    }
    double ans = 0;
    for (auto tmp : _coeff) {
      ans = ans * x + tmp;
    }
    // return ans;
    return fabs(ans - 0.0) < EPSILON ? 0.0 : ans;
  }
  double gradientAt(double x) {
    double ans = 0;
    for (auto tmp : _gradient) {
      ans = ans * x + tmp;
    }
    return ans;
  }

  // Overload operation
  const double& operator[](int i) const { return _coeff[i]; }
  Poly operator+(Poly& b) {
    vector<double> d1 = _coeff, d2 = b.getCoef(),
                   r(fmax(size(), b.size()), 0.0);
    reverse(d1.begin(), d1.end());
    reverse(d2.begin(), d2.end());
    for (int i = 0; i < d1.size(); i++) r[i] += d1[i];
    for (int i = 0; i < d2.size(); i++) r[i] += d2[i];
    reverse(r.begin(), r.end());
    return Poly(r);
  }
  Poly& operator+=(Poly& b) {
    vector<double> d1 = _coeff, d2 = b.getCoef(),
                   r(fmax(size(), b.size()), 0.0);
    reverse(d1.begin(), d1.end());
    reverse(d2.begin(), d2.end());
    for (int i = 0; i < d1.size(); i++) r[i] += d1[i];
    for (int i = 0; i < d2.size(); i++) r[i] += d2[i];
    reverse(r.begin(), r.end());
    _coeff = r;
    _gradCalc();
    _cleanUp();
    return *this;
  }
  Poly operator-(Poly& b) {
    vector<double> d1 = _coeff, d2 = b.getCoef(),
                   r(fmax(size(), b.size()), 0.0);
    reverse(d1.begin(), d1.end());
    reverse(d2.begin(), d2.end());
    for (int i = 0; i < d1.size(); i++) r[i] += d1[i];
    for (int i = 0; i < d2.size(); i++) r[i] -= d2[i];
    reverse(r.begin(), r.end());
    Poly ret = Poly(r);
    // return Poly(r);
    return ret;
  }
  Poly& operator-=(Poly& b) {
    vector<double> d1 = _coeff, d2 = b.getCoef(),
                   r(fmax(size(), b.size()), 0.0);
    reverse(d1.begin(), d1.end());
    reverse(d2.begin(), d2.end());
    for (int i = 0; i < d1.size(); i++) r[i] += d1[i];
    for (int i = 0; i < d2.size(); i++) r[i] -= d2[i];
    reverse(r.begin(), r.end());
    _coeff = r;
    _gradCalc();
    _cleanUp();
    return *this;
  }
  Poly operator*(Poly& b) {
    int N = b.size() + size() - 1, N1 = size(), N2 = b.size();
    vector<double> r(N, 0.0);
    int idx1, idx2;
    for (int i = 0; i < N1; i++) {
      idx1 = N1 - i - 1;
      for (int j = 0; j < N2; j++) {
        idx2 = N2 - j - 1;
        r[N - 1 - idx1 - idx2] += _coeff[i] * b[j];
      }
    }
    return Poly(r);
  }
  Poly& operator*=(Poly& b) {
    int N = b.size() + size() - 1, N1 = size(), N2 = b.size();
    vector<double> r(N, 0.0);
    int idx1, idx2;
    for (int i = 0; i < N1; i++) {
      idx1 = N1 - i - 1;
      for (int j = 0; j < N2; j++) {
        idx2 = N2 - j - 1;
        r[N - 1 - idx1 - idx2] += _coeff[i] * b[j];
      }
    }
    _coeff = r;
    _gradCalc();
    _cleanUp();
    return *this;
  }
  Poly operator/(Poly& b) {
    vector<double> tmp(b.size(), 0.0);
    Poly q = Poly(tmp), r = *this, t, tmpp;
    while (!r.isZero() && r.deg() >= b.deg()) {
      t = r.leadDiv(b);
      q += t;
      tmpp = t * b;
      r -= tmpp;
    }
    return q;
  }
  Poly& operator/=(Poly& b) {
    vector<double> tmp(b.size(), 0.0);
    Poly q = Poly(tmp), r = *this, t, tmpp;
    while (!r.isZero() && r.deg() >= b.deg()) {
      t = r.leadDiv(b);
      q += t;
      tmpp = t * b;
      r -= tmpp;
    }
    _coeff = q.getCoef();
    _gradCalc();
    _cleanUp();
    return (*this);
  }
  Poly operator*(const double t) {
    vector<double> d = _coeff;
    for (int i = 0; i < d.size(); i++) d[i] *= t;
    return Poly(d);
  }
};

static std::ostream& operator<<(std::ostream& out, const Poly& u) {
  int N = u.size();
  for (int i = 0; i < N; i++) {
    if (u[i] == 0) continue;
    char sign = '+';
    if (u[i] < 0) sign = '-';
    out << " " << sign << " ";
    out << fabs(u[i]) << 'x' << '^' << N - i - 1;
  }
  return out;
}

#endif
