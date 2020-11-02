#ifndef BUDAN_H
#define BUDAN_H

#include "util.h"

vector<tuple<double, double>> budanIsoroot(Poly& p) {
    vector<tuple<double, double>> ret;
    deque<Boundry> b;
    double tmp = bound(p), mid, root;
    int midchange;
    Boundry tmpb;

    b.push_back(
        Boundry{-tmp, signChangeNum(p, -tmp), tmp, signChangeNum(p, tmp)});

    while (b.size() > 0) {
      tmpb = b[0];
      b.pop_front();
      if (DEBUG_BUDAN)
        cout << "Search Boundary: " << tmpb.left << "\t | \t  " << tmpb.right
             << endl;
      if (DEBUG_BUDAN)
        cout << tmpb.lchange << "\t | \t  " << tmpb.rchange << endl;
      mid = (tmpb.left + tmpb.right) / 2;
      midchange = signChangeNum(p, mid);
      if (DEBUG_BUDAN) cout << "Mid Change:  " << midchange << endl;

      // left side
      if (mid - tmpb.left < MINRANGE && tmpb.lchange - midchange > 0) {
        root = rootInBound(p, tmpb.left, mid);
        if (root != NOTFOUND)
          ret.push_back(tuple<double, double>(tmpb.left, mid));
      } else if (tmpb.lchange - midchange > 0) {
        b.push_back(Boundry{tmpb.left, tmpb.lchange, mid, midchange});
      }

      // right side
      if (tmpb.right - mid < MINRANGE && midchange - tmpb.rchange > 0) {
        root = rootInBound(p, mid, tmpb.right);
        if (root != NOTFOUND)
          ret.push_back(tuple<double, double>(mid, tmpb.right));
      } else if (midchange - tmpb.rchange > 0) {
        b.push_back(Boundry{mid, midchange, tmpb.right, tmpb.rchange});
      }
    }

    return ret;
  }

vector<double> budanSolve(Poly& p) {
  if (p.deg() == 1) return {-p[1] / p[0]};
  vector<Poly> plist = squareFreeDecompo(p);
  vector<tuple<double, double>> ranges, tmprange;
  vector<double> roots;
  double tmp;

  if (DEBUG_BUDAN) {
    cout << "DEBUG BUDAN SQUREFREE DECOMP: " << endl;
    for (auto x : plist) cout << x << endl;
  }

  for (auto p : plist) {
    if (p.deg() <= 0) continue;
    tmprange = budanIsoroot(p);
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_BUDAN)
      cout << "Finding root in range: " << get<0>(range) << " to "
           << get<1>(range) << endl;
    tmp = rootInBound(p, get<0>(range), get<1>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }

  return roots;
}

#endif
