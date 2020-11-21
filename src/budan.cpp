#include "budan.h"

#include <iostream>

using std::cout;
using std::deque;
using std::endl;

vector<tuple<double, double>> budanIsoroot(Poly& p) {
  vector<tuple<double, double>> ret;
  deque<Boundry> b;
  double tmp = upperBound(p), mid, root;
  int midchange;
  Boundry tmpb;

  b.push_back(
      Boundry{-tmp, signChangeNum(p, -tmp), tmp, signChangeNum(p, tmp)});

  while (b.size() > 0) {
    tmpb = b[0];
    b.pop_front();
    mid = (tmpb.left + tmpb.right) / 2;
    midchange = signChangeNum(p, mid);

    // left side
    if (mid - tmpb.left < MINRANGE && tmpb.lchange - midchange > 0) {
      ret.push_back(tuple<double, double>(tmpb.left, mid));
    } else if (tmpb.lchange - midchange > 0) {
      b.push_back(Boundry{tmpb.left, tmpb.lchange, mid, midchange});
    }

    // right side
    if (tmpb.right - mid < MINRANGE && midchange - tmpb.rchange > 0) {
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

  for (auto p : plist) {
    if (DEBUG_BUDAN) cout << p << endl;
    if (p.deg() <= 0) continue;
    tmprange = budanIsoroot(p);
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_BUDAN)
      std::cout << "DEBUG_BUDAN: Searching root in range: "
                << std::get<0>(range) << " to " << std::get<1>(range)
                << std::endl;
    tmp = rootInBound(p, std::get<0>(range), std::get<1>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }

  return roots;
}
