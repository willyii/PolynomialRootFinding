#include "vincent.h"

vector<tuple<double, double>> Vincent::isoRoot(Poly& p) {
  vector<tuple<double, double>> ret;
  double maxOfRoot = bound(p);

  // Init L
  deque<tuple<Poly, Poly>> L;
  vector<double> tmpCoef = {1, 0};
  Poly tmpPoly = Poly(tmpCoef);
  L.emplace_back(tuple<Poly, Poly>(p, tmpPoly));
  tmpCoef[0] = -1;
  tmpPoly = Poly(tmpCoef);
  L.emplace_back(tuple<Poly, Poly>(timeToP(p, -1), tmpPoly));

  while (L.size() > 0) {
    auto frontTuple = L.front();
    L.pop_front();
    Poly A = get<0>(frontTuple), M = get<1>(frontTuple);
    int v = signChangeNum(A, 0);
    if (v == 0) continue;
    if (v == 1)
      ret.emplace_back(tuple<double, double>(
          M.valueAt(0),
          fmax(-maxOfRoot, fmin(maxOfRoot, M.valueAt(INFINITY)))));

    // TODO: Optimize b
    int b = 1;
    Poly B = addToP(A, b);
    int w = v - signChangeNum(B, 0);
    if (B.valueAt(0) == 0) {
      ret.emplace_back(tuple<double, double>(M.valueAt(b), M.valueAt(b)));
      tmpCoef = {1, 0};
      tmpPoly = Poly(tmpCoef);
      B /= tmpPoly;
    }
    L.emplace_back(tuple<Poly, Poly>(B, addToP(M, b)));
    if (w == 0) continue;
    if (w == 1) ret.emplace_back(M.valueAt(0), M.valueAt(b));
    if (w > 1) {
      L.emplace_back(
          tuple<Poly, Poly>(inverseTimesB(A, b), inverseTimesB(M, b)));
    }
  }
  return ret;
}

vector<double> Vincent::solve(Poly& p) {
  if (p.deg() == 1) return {-p[1] / p[0]};
  vector<Poly> plist = squareFreeDecompo(p);
  vector<tuple<double, double>> ranges, tmprange;
  set<double> roots;
  double tmp;

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
    if (tmp != NOTFOUND) roots.emplace(tmp);
  }

  vector<double> ans;
  ans.assign(roots.begin(), roots.end());
  return ans;
}

Poly Vincent::inverseTimesB(Poly& p, int b) {
  int N = p.deg();
  vector<double> tmpAdd1 = {1, 1};
  vector<double> tmpB = {double(b)};
  Poly util = Poly(tmpB);
  Poly pB = util.pow(N);
  Poly addOne = Poly(tmpAdd1);
  vector<double> init = {0};

  Poly ans = pB * p[0];
  pB /= util;

  for (size_t i = 1; i < N + 1; i++) {
    Poly tmp = pB * addOne;
    tmp = tmp * p[i];
    ans += tmp;
    pB /= util;
    addOne *= addOne;
  }
  ans /= addOne;
  return ans;
}