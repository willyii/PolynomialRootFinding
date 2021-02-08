// ----------------------------------------------------------------------------
//
// FILENAME: interval.h
//
// DESCRIPTION:
//    This file define a class represent for the interval that used in interval
//    arithmetic.
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_INTERVAL_H_
#define POLY_INTERVAL_H_

#include <algorithm>
#include <fstream>
#include <limits>

#include "param.h"

/* TODO : Use library to do interval */
static const double kROUNDOFF = 1e-13;

struct interval {
  double left;
  double right;

  interval() : left(-kROUNDOFF), right(kROUNDOFF){};
  interval(double num) : left(num - kROUNDOFF), right(num + kROUNDOFF){};
  interval(double num1, double num2) : left(num1), right(num2){};

  double value() const {
    if (left < 0 && right > 0)
      return 0.0;
    return (left + right) / 2;
    // return left;
  }

  inline interval &operator=(interval a) {
    left = a.left;
    right = a.right;
    return *this;
  }

  inline interval &operator=(double a) {
    left = a - kROUNDOFF;
    right = a + kROUNDOFF;
    return *this;
  }

  inline interval &operator+=(interval a) {
    left += a.left;
    right += a.right;
    return *this;
  }

  inline interval &operator+=(double a) {
    left += a;
    right += a;
    return *this;
  }

  inline interval &operator-=(interval a) {
    left -= a.right;
    right -= a.left;
    return *this;
  }

  inline interval &operator-=(double a) {
    left -= a;
    right -= a;
    return *this;
  }

  inline interval &operator*=(interval a) {
    double num1(left * a.left);
    double num2(left * a.right);
    double num3(right * a.left);
    double num4(right * a.right);

    left = std::min({num1, num2, num3, num4});
    right = std::max({num1, num2, num3, num4});

    return *this;
  }

  inline interval &operator*=(double a) {
    left *= a;
    right *= a;
    if (a < 0)
      std::swap(left, right);
    return *this;
  }

  inline interval &operator/=(interval a) {
    // if 0 in a
    if (a.left <= 0 && a.right >= 0) {
      left = (left > 0.0) ? std::numeric_limits<double>::infinity()
                          : -std::numeric_limits<double>::infinity();
      right = (right > 0.0) ? std::numeric_limits<double>::infinity()
                            : -std::numeric_limits<double>::infinity();
      return *this;
    }
    left /= a.right;
    right /= a.left;
    return *this;
  }

  inline interval &operator/=(double a) {
    left /= a;
    right /= a;
    if (a < 0) {
      std::swap(left, right);
    }
    return *this;
  }
};

interval operator+(const interval &interval1, const interval &interval2) {
  return interval(interval1) += interval2;
}

interval operator+(const interval &interval1, double num) {
  return interval(interval1) += num;
}

interval operator+(double num, const interval &interval1) {
  return interval(interval1) += num;
}

interval operator-(const interval &interval1, const interval &interval2) {
  return interval(interval1) -= interval2;
}

interval operator-(const interval &interval1, double num) {
  return interval(interval1) -= num;
}

interval operator-(double num, const interval &interval1) {
  return (interval(interval1) *= -1.0) += num;
}

interval operator*(const interval &interval1, const interval &interval2) {
  return interval(interval1) *= interval2;
}

interval operator*(const interval &interval1, double num) {
  return interval(interval1) *= num;
}

interval operator*(double num, const interval &interval1) {
  return interval(interval1) *= num;
}

interval operator/(const interval &interval1, const interval &interval2) {
  return interval(interval1) /= interval2;
}

interval operator/(const interval &interval1, double num) {
  return interval(interval1) /= num;
}

interval operator/(double num, const interval &interval1) {
  interval tmp(num);
  return interval(tmp) /= interval1;
}
#endif // POLY_INTERVAL_H_
