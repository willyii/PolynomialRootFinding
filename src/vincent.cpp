#include "vincent.h"

vector<tuple<double, double>> Vincent::isoRoot(Poly& p) {
    return {}
}


vector<double> Vincent::solve(Poly& p) {
  if (p.deg() == 1) return {-p[1] / p[0]};
  vector<Poly> plist = squareFreeDecompo(p);
  vector<tuple<double, double>> ranges, tmprange;
  vector<double> roots;
  double tmp;

  if (DEBUG_VINCENT) {
    cout << "DEBUG Vincent SQUREFREE DECOMP: " << endl;
    for (auto x : plist) cout << x << endl;
  }

  for (auto p : plist) {
    if (p.deg() <= 0) continue;
    tmprange = isoRoot(p);
    ranges.insert(ranges.end(), tmprange.begin(), tmprange.end());
  }

  for (auto range : ranges) {
    if (DEBUG_VINCENT)
      cout << "Finding root in range: " << get<0>(range) << " to "
           << get<1>(range) << endl;
    tmp = rootInBound(p, get<0>(range), get<1>(range));
    if (tmp != NOTFOUND) roots.emplace_back(tmp);
  }

  return roots;
}