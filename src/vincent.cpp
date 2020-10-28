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
      ret.emplace_back(
          tuple<double, double>(M.valueAt(0),fmax(-maxOfRoot, fmin(maxOfRoot, M.valueAt(INFINITY))))
          );
    
    // TODO: Optimize b
    int b = 1;
    Poly B = addToP(A, b);
    int w = v - signChangeNum(B, 0);
    if(B.valueAt(0) == 0){
      ret.emplace_back(tuple<double, double>(M.valueAt(b), M.valueAt(b)));
      tmpCoef = {1, 0};
      tmpPoly = Poly(tmpCoef);
      B /= tmpPoly;
    }

    L.emplace_back(tuple<Poly, Poly>(B, addToP(M, b)));
    if(w == 0) continue;
    if(w==1) ret.emplace_back(M.valueAt(0), M.valueAt(b));
    if(w > 1) continue; // TODO: util and reload ^ operator

  }

  return {};
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

Poly Vincent::inverseTimesB(Poly& p, int b){

  

  return Poly();
}