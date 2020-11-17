#include "vincent.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <tuple>

#include "util.h"

using std::deque;
using std::get;
using std::make_tuple;
using std::reverse;
using std::swap;

vector<tuple<double, double>> vincentIsoroot(Poly& p, bool negative) {
  vector<tuple<double, double>> rootlist;
  int s;
  double upper = upperBound(p), lower = lowerBound(p);
  s = signChangeNum(p, 0);
  double a, b, c, d, a0 = 16;
  Poly q;
  if (s == 0) return {};
  if (s == 1) return {make_tuple(0.0, upper)};
  deque<Interval> intervals = {Interval{1, 0, 0, 1, p, s}};
  // Helper poly x;
  vector<double> tmp = {1, 0};
  Poly tmpPoly = Poly(tmp);
  while (intervals.size() > 0) {
    // pop from front
    auto current = intervals[0];
    intervals.pop_front();
    a = current.a;
    b = current.b;
    c = current.c;
    d = current.d;
    q = current.p;
    s = current.s;
    lower = lowerBound(q);
    if (lower > a0) {
      q = timeToP(q, lower);
      a *= lower;
      c *= lower;
      lower = 1;
    }
    if (lower >= 1) {
      q = addToP(q, lower);
      b += lower * a;
      c *= lower;
    }
    if (q.valueAt(0) == 0.0) {
      if (DEBUG_VINCENT)
        std::cout << "DEBUG_VINCENT: b/d = " << b / d << std::endl;
      rootlist.emplace_back(make_tuple(b / d, b / d));
      q /= tmpPoly;
      s = signChangeNum(q, 0);
      if (s == 0) continue;
    }

    Poly q1 = addToP(q, 1);
    double a1 = a, b1 = a + b, c1 = c, d1 = c + d;
    int r = 0;

    if (q1.valueAt(0) == 0.0) {
      rootlist.emplace_back(make_tuple(b1 / d1, b1 / d1));
      q1 /= tmpPoly;
      r = 1;
    }
    int s1 = signChangeNum(q1, 0);
    int s2 = s - s1 - r;
    double a2 = b, b2 = a + b, c2 = d, d2 = c + d;
    Poly q2;
    if (s2 > 1) {
      q2 = reverse(q);
      q2 = addToP(q2, 1);
      if (q2.valueAt(0) == 0) q2 /= tmpPoly;
      s2 = signChangeNum(q2, 0);
    }

    if (s1 < s2) {
      swap(a1, a2);
      swap(b1, b2);
      swap(c1, c2);
      swap(d1, d2);
      swap(q1, q2);
      swap(s1, s2);
    }

    if (s1 == 0) continue;
    if (s1 == 1) {
      // if (DEBUG_VINCENT)
      //  std::cout << "DEBUG_VINCENT: a1/c1 = " << a1 / c1
      //            << ", b1 / d1 = " << b1 / d1 << std::endl;
      // if (DEBUG_VINCENT)
      //  std::cout << "DEBUG_VINCENT: a1 = " << a1 << ", b1 = " << b1
      //            << ", c1 = " << c1 << ", d1 = " << d1 << std::endl;
      double start = fmin(fmin(a1 / c1, upper), fmin(b1 / d1, upper));
      double end = fmax(fmin(a1 / c1, upper), fmin(b1 / d1, upper));
      rootlist.emplace_back(make_tuple(start, end));
      // if (fabs(c1 - 0) < EPSILON)
      //  rootlist.emplace_back(make_tuple(b1 / d1, upper));
      // else
      //  rootlist.emplace_back(make_tuple(a1 / c1, b1 / d1));
    } else {
      intervals.emplace_back(Interval{a1, b1, c1, d1, q1, s1});
    }
    if (s2 == 0) continue;
    if (s2 == 1) {
      // if (DEBUG_VINCENT)
      //  std::cout << "DEBUG_VINCENT: a2/c2 = " << a2 / c2
      //            << ", b2 / d2 = " << b2 / d2 << std::endl;
      double start = fmin(b2 / d2, a2 / c2);
      double end = fmax(b2 / d2, a2 / c2);
      rootlist.emplace_back(make_tuple(start, end));
    } else {
      intervals.emplace_back(Interval{a2, b2, c2, d2, q2, s2});
    }
  }

  // if (negative == true) {
  //  for (size_t i = 0; i < rootlist.size(); i++) {
  //    std::cout << "Last Bug BE: " << get<0>(rootlist[i]) << " and "
  //              << get<1>(rootlist[i]) << std::endl;
  //    double tmp = get<0>(rootlist[i]);
  //    get<0>(rootlist[i]) = -get<1>(rootlist[i]);
  //    get<1>(rootlist[i]) = -tmp;
  //    // rootlist[i] = std::make_tuple(-get<1>(rootlist[i]),
  //    // -get<0>(rootlist[i]));
  //    std::cout << "Last Bug AF: " << get<0>(rootlist[i]) << " and "
  //              << get<1>(rootlist[i]) << std::endl;
  //  }
  //}

  return rootlist;
}

vector<double> vincentSolve(Poly& p) {
  if (p.deg() == 1) return {-p[1] / p[0]};
  vector<Poly> plist = squareFreeDecompo(p);
  vector<tuple<double, double>> ranges, tmprange;
  vector<double> roots;
  double tmp;

  // Postive roots
  for (auto poly : plist) {
    if (poly.deg() <= 0) continue;
    tmprange = vincentIsoroot(poly, false);
    for (size_t i = 0; i < tmprange.size(); i++) {
      // if (DEBUG_VINCENT)
      //  std::cout << "Before Range: from " << get<0>(tmprange[i]) << " to "
      //            << get<1>(tmprange[i]) << std::endl;
      refineRange(poly, tmprange[i]);
      // if (DEBUG_VINCENT)
      //  std::cout << "After Range: from " << get<0>(tmprange[i]) << " to "
      //            << get<1>(tmprange[i]) << std::endl;
    }
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_VINCENT)
      std::cout << "DEBUG VINCENT: Searching root in range from "
                << std::get<0>(range) << " to " << std::get<1>(range)
                << std::endl;
    tmp = rootInBound(p, std::get<0>(range), std::get<1>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }

  // Negative roots
  ranges = {};
  for (auto poly : plist) {
    if (poly.deg() <= 0) continue;
    Poly tmpp = timeToP(poly, -1);
    tmprange = vincentIsoroot(tmpp, false);
    for (size_t i = 0; i < tmprange.size(); i++) {
      // if (DEBUG_VINCENT)
      //  std::cout << "Before Range: from " << -get<1>(tmprange[i]) << " to "
      //            << -get<0>(tmprange[i]) << std::endl;
      refineRange(tmpp, tmprange[i]);
      // if (DEBUG_VINCENT)
      //  std::cout << "After Range: from " << -get<1>(tmprange[i]) << " to "
      //            << -get<0>(tmprange[i]) << std::endl;
    }
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_VINCENT)
      std::cout << "DEBUG VINCENT: Searching root in range from "
                << -std::get<1>(range) << " to " << -std::get<0>(range)
                << std::endl;
    tmp = rootInBound(p, -std::get<1>(range), -std::get<0>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }
  return roots;
}

Poly reverse(Poly& p) {
  vector<double> coef = p.getCoef();
  reverse(coef.begin(), coef.end());
  Poly ans = Poly(coef);
  return ans;
}

// Refine the range, make the range smaller enough to avoid edge case
void refineRange(Poly& p, tuple<double, double>& range) {
  int lchange = signChangeNum(p, get<0>(range));
  int rchange = signChangeNum(p, get<1>(range));
  while ((get<1>(range) - get<0>(range)) > MINRANGE) {
    double mid = (get<0>(range) + get<1>(range)) / 2;
    int midchange = signChangeNum(p, mid);
    if ((lchange - midchange) % 2 == 0) {
      get<0>(range) = mid;
      lchange = midchange;
    } else {
      get<1>(range) = mid;
      rchange = midchange;
    }
  }
  return;
}
