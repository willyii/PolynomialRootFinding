#ifndef COEF_H
#define COEF_H

#include <math.h>

#include <vector>

#include "param.h"

using std::vector;

struct Coef {
  vector<double> data;

  Coef() {}

  Coef(const vector<double>& d) { data = d; }

  Coef(int n) { data = vector<double>(n, 0.0); }

  int size() const { return data.size(); }

  double lc() const {
    for (int i = 0; i < data.size(); i++)
      if (abs(data[i] - 0) > EPSILON) return data[i];
    return 0;
  }

  int deg() const {
    int N = data.size();
    for (int i = 0; i < N; i++) {
      if (abs(data[i] - 0) > EPSILON) {
        return N - i - 1;
      }
    }
    return -1;
  }

  void fillZero() {
    for (int i = 0; i < data.size(); i++) data[i] = 0;
  }

  Coef operator+(const Coef& d) {
    Coef r = Coef(d.size());
    for (int i = 0; i < data.size(); i++) r[i] = data[i] + d[i];
    return r;
  }

  Coef& operator+=(const Coef& d) {
    for (int i = 0; i < data.size(); i++) data[i] += d[i];
    return *this;
  }

  Coef operator-(const Coef& d) {
    Coef r = Coef(d.size());
    for (int i = 0; i < data.size(); i++) r[i] = data[i] - d[i];
    return r;
  }

  Coef& operator-=(const Coef& d) {
    for (int i = 0; i < data.size(); i++) data[i] -= d[i];
    return *this;
  }

  Coef operator*(const Coef& d) {
    Coef r = Coef(data.size());
    int idx1, idx2, N = data.size();
    for (int i = 0; i < N; i++) {
      idx1 = N - i - 1;
      for (int j = 0; j < N; j++) {
        idx2 = N - j - 1;
        if (idx1 + idx2 >= N) continue;
        r[N - 1 - idx1 - idx2] += data[i] * d[j];
      }
    }
    return r;
  }

  Coef& operator*=(const Coef& d) {
    vector<double> r(data.size(), 0);
    int idx1, idx2, N = data.size();
    for (int i = 0; i < N; i++) {
      idx1 = N - i - 1;
      for (int j = 0; j < N; j++) {
        idx2 = N - j - 1;
        if (idx1 + idx2 >= N) continue;
        r[N - 1 - idx1 - idx2] += data[i] * d[j];
      }
    }
    data = r;
    return *this;
  }

  Coef operator/(const Coef& d) {
    Coef q = Coef(d.size()), r = *this, t;
    while (!r.isZero() && r.deg() >= d.deg()) {
      t = r.leadDiv(d);
      q += t;
      r -= t * d;
    }
    return q;
  }

  Coef& operator/=(const Coef& d) {
    Coef q = Coef(d.size()), r = *this, t;
    while (!r.isZero() && r.deg() >= d.deg()) {
      t = r.leadDiv(d);
      q += t;
      r -= t * d;
    }
    (*this) = q;
    return (*this);
  }

  const double& operator[](int i) const { return data[i]; }

  double& operator[](int i) { return data[i]; }

  /*----------------------------utils----------------------------*/
  bool isZero() {
    for (int i = 0; i < data.size(); i++)
      if (abs(data[i] - 0) > EPSILON) return false;
    return true;
  }

  Coef leadDiv(const Coef& c2) {
    double div = this->lc() / c2.lc();
    int degree = this->deg() - c2.deg();
    Coef t = Coef(data.size());
    t[data.size() - degree - 1] = div;
    return t;
  }
};

static std::ostream& operator<<(std::ostream& out, const Coef& u) {
  for (int i = 0; i < u.size(); i++) {
    if (i) out << ' ';
    out << u[i];
  }
  return out;
}

typedef Coef coef;

#endif